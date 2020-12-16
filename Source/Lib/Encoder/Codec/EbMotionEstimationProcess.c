/*
* Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#include <stdlib.h>

#include "EbEncHandle.h"
#include "EbUtility.h"
#include "EbPictureControlSet.h"
#include "EbPictureDecisionResults.h"
#include "EbMotionEstimationProcess.h"
#include "EbMotionEstimationResults.h"
#include "EbReferenceObject.h"
#include "EbMotionEstimation.h"
#include "EbLambdaRateTables.h"
#include "EbComputeSAD.h"
#ifdef ARCH_X86_64
#include <emmintrin.h>
#endif
#include "EbTemporalFiltering.h"
#include "EbGlobalMotionEstimation.h"

#include "EbResize.h"
#include "EbPictureDemuxResults.h"
#include "EbRateControlTasks.h"
#include "firstpass.h"
#include "EbInitialRateControlProcess.h"
/* --32x32-
|00||01|
|02||03|
--------*/
/* ------16x16-----
|00||01||04||05|
|02||03||06||07|
|08||09||12||13|
|10||11||14||15|
----------------*/
/* ------8x8----------------------------
|00||01||04||05|     |16||17||20||21|
|02||03||06||07|     |18||19||22||23|
|08||09||12||13|     |24||25||28||29|
|10||11||14||15|     |26||27||30||31|

|32||33||36||37|     |48||49||52||53|
|34||35||38||39|     |50||51||54||55|
|40||41||44||45|     |56||57||60||61|
|42||43||46||47|     |58||59||62||63|
-------------------------------------*/
EbErrorType check_00_center(PictureParentControlSet *pcs_ptr, EbPictureBufferDesc *ref_pic_ptr,
                            MeContext *context_ptr, uint32_t sb_origin_x, uint32_t sb_origin_y,
                            uint32_t sb_width, uint32_t sb_height, int16_t *x_search_center,
                            int16_t *y_search_center);

/************************************************
 * Set ME/HME Params from Config
 ************************************************/
void *set_me_hme_params_from_config(SequenceControlSet *scs_ptr, MeContext *me_context_ptr) {
    uint16_t hme_region_index = 0;

    me_context_ptr->search_area_width  = (uint8_t)scs_ptr->static_config.search_area_width;
    me_context_ptr->search_area_height = (uint8_t)scs_ptr->static_config.search_area_height;

    me_context_ptr->number_hme_search_region_in_width =
        (uint16_t)scs_ptr->static_config.number_hme_search_region_in_width;
    me_context_ptr->number_hme_search_region_in_height =
        (uint16_t)scs_ptr->static_config.number_hme_search_region_in_height;

    me_context_ptr->hme_level0_total_search_area_width =
        (uint16_t)scs_ptr->static_config.hme_level0_total_search_area_width;
    me_context_ptr->hme_level0_total_search_area_height =
        (uint16_t)scs_ptr->static_config.hme_level0_total_search_area_height;

    for (hme_region_index = 0; hme_region_index < me_context_ptr->number_hme_search_region_in_width;
         ++hme_region_index) {
        me_context_ptr->hme_level0_search_area_in_width_array[hme_region_index] =
            (uint16_t)
                scs_ptr->static_config.hme_level0_search_area_in_width_array[hme_region_index];
        me_context_ptr->hme_level1_search_area_in_width_array[hme_region_index] =
            (uint16_t)
                scs_ptr->static_config.hme_level1_search_area_in_width_array[hme_region_index];
        me_context_ptr->hme_level2_search_area_in_width_array[hme_region_index] =
            (uint16_t)
                scs_ptr->static_config.hme_level2_search_area_in_width_array[hme_region_index];
    }

    for (hme_region_index = 0;
         hme_region_index < me_context_ptr->number_hme_search_region_in_height;
         ++hme_region_index) {
        me_context_ptr->hme_level0_search_area_in_height_array[hme_region_index] =
            (uint16_t)
                scs_ptr->static_config.hme_level0_search_area_in_height_array[hme_region_index];
        me_context_ptr->hme_level1_search_area_in_height_array[hme_region_index] =
            (uint16_t)
                scs_ptr->static_config.hme_level1_search_area_in_height_array[hme_region_index];
        me_context_ptr->hme_level2_search_area_in_height_array[hme_region_index] =
            (uint16_t)
                scs_ptr->static_config.hme_level2_search_area_in_height_array[hme_region_index];
    }

    return NULL;
}

/************************************************
 * Set ME/HME Params
 ************************************************/
void *set_me_hme_params_oq(MeContext *me_context_ptr, PictureParentControlSet *pcs_ptr,
                           SequenceControlSet *scs_ptr, EbInputResolution input_resolution) {
    UNUSED(scs_ptr);
    // HME/ME default settings
    me_context_ptr->number_hme_search_region_in_width  = 2;
    me_context_ptr->number_hme_search_region_in_height = 2;

    // Set the minimum ME search area
    if (pcs_ptr->sc_content_detected)
            if (pcs_ptr->enc_mode <= ENC_M3) {
                me_context_ptr->search_area_width = me_context_ptr->search_area_height = 175;
                me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 750;
            }
            else if (pcs_ptr->enc_mode <= ENC_M5) {
                me_context_ptr->search_area_width = me_context_ptr->search_area_height = 125;
                me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 500;
            }
            else {
                if (use_output_stat(scs_ptr) || (scs_ptr->lap_enabled && !pcs_ptr->first_pass_done)) {
                    me_context_ptr->search_area_width = me_context_ptr->search_area_height = 37;
                    me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 175;
                }
                else {
                me_context_ptr->search_area_width = me_context_ptr->search_area_height = 75;
                me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 350;
        }
            }
    else if (pcs_ptr->enc_mode <= ENC_M0) {
    me_context_ptr->search_area_width = me_context_ptr->search_area_height = 64;
    me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 256;
    }
    else if (pcs_ptr->enc_mode <= ENC_M1) {
        me_context_ptr->search_area_width = me_context_ptr->search_area_height = 64;
        me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 192;
    }
    else if (pcs_ptr->enc_mode <= ENC_M2) {
        me_context_ptr->search_area_width = me_context_ptr->search_area_height = 64;
        me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 128;
    }
    else if (pcs_ptr->enc_mode <= ENC_M5) {
        if (use_output_stat(scs_ptr)) {
            me_context_ptr->search_area_width = me_context_ptr->search_area_height = 8;
            me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 8;
        }
        else {
        me_context_ptr->search_area_width = me_context_ptr->search_area_height = 16;
        me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 64;
        }
    }
    else {
        if (use_output_stat(scs_ptr) || (scs_ptr->lap_enabled && !pcs_ptr->first_pass_done)) {
            me_context_ptr->search_area_width = me_context_ptr->search_area_height = 8;
            me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 8;
        }
        else {
            me_context_ptr->search_area_width = me_context_ptr->search_area_height = 16;
            me_context_ptr->max_me_search_width = 64;
            me_context_ptr->max_me_search_height = 32;
        }
    }
        if (pcs_ptr->enc_mode <= ENC_M1) {
            me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = input_resolution <= INPUT_SIZE_1080p_RANGE ? 120 : 240;
            me_context_ptr->hme_level0_max_total_search_area_width = me_context_ptr->hme_level0_max_total_search_area_height = 480;
        }
        else {
            me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = 32;
            me_context_ptr->hme_level0_max_total_search_area_width = me_context_ptr->hme_level0_max_total_search_area_height = 164;
        }
    if (!pcs_ptr->sc_content_detected)
        if (use_output_stat(scs_ptr) || (scs_ptr->lap_enabled && !pcs_ptr->first_pass_done)) {
            me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height  = me_context_ptr->hme_level0_total_search_area_width/2;
            me_context_ptr->hme_level0_max_total_search_area_width = me_context_ptr->hme_level0_max_total_search_area_height =   me_context_ptr->hme_level0_max_total_search_area_width/2;
        }
    me_context_ptr->hme_level0_max_search_area_in_width_array[0] =
        me_context_ptr->hme_level0_max_search_area_in_width_array[1] =
        me_context_ptr->hme_level0_max_total_search_area_width / me_context_ptr->number_hme_search_region_in_width;
    me_context_ptr->hme_level0_max_search_area_in_height_array[0] =
        me_context_ptr->hme_level0_max_search_area_in_height_array[1] =
        me_context_ptr->hme_level0_max_total_search_area_height / me_context_ptr->number_hme_search_region_in_height;
    me_context_ptr->hme_level0_search_area_in_width_array[0] =
        me_context_ptr->hme_level0_search_area_in_width_array[1] =
        me_context_ptr->hme_level0_total_search_area_width / me_context_ptr->number_hme_search_region_in_width;
    me_context_ptr->hme_level0_search_area_in_height_array[0] =
        me_context_ptr->hme_level0_search_area_in_height_array[1] =
        me_context_ptr->hme_level0_total_search_area_height / me_context_ptr->number_hme_search_region_in_height;
    if (pcs_ptr->enc_mode <= ENC_M4) {
        me_context_ptr->hme_level1_search_area_in_width_array[0] =
            me_context_ptr->hme_level1_search_area_in_width_array[1] =
            me_context_ptr->hme_level1_search_area_in_height_array[0] =
            me_context_ptr->hme_level1_search_area_in_height_array[1] = 16;
    }
    else {
        me_context_ptr->hme_level1_search_area_in_width_array[0] =
            me_context_ptr->hme_level1_search_area_in_width_array[1] = 8;
        me_context_ptr->hme_level1_search_area_in_height_array[0] =
            me_context_ptr->hme_level1_search_area_in_height_array[1] = 3;
    }
    if (pcs_ptr->enc_mode <= ENC_M4) {
        me_context_ptr->hme_level2_search_area_in_width_array[0] =
            me_context_ptr->hme_level2_search_area_in_width_array[1] =
            me_context_ptr->hme_level2_search_area_in_height_array[0] =
            me_context_ptr->hme_level2_search_area_in_height_array[1] = 16;
    }
    else {
        me_context_ptr->hme_level2_search_area_in_width_array[0] =
            me_context_ptr->hme_level2_search_area_in_width_array[1] = 8;

        me_context_ptr->hme_level2_search_area_in_height_array[0] =
            me_context_ptr->hme_level2_search_area_in_height_array[1] = 3;
    }
    if (!pcs_ptr->sc_content_detected)
        if (use_output_stat(scs_ptr) || (scs_ptr->lap_enabled && !pcs_ptr->first_pass_done)) {
            me_context_ptr->hme_level1_search_area_in_width_array[0] =
                me_context_ptr->hme_level1_search_area_in_width_array[1] =
                me_context_ptr->hme_level1_search_area_in_height_array[0] =
                me_context_ptr->hme_level1_search_area_in_height_array[1] = 16/2;

            me_context_ptr->hme_level2_search_area_in_width_array[0] =
                me_context_ptr->hme_level2_search_area_in_width_array[1] =
                me_context_ptr->hme_level2_search_area_in_height_array[0] =
                me_context_ptr->hme_level2_search_area_in_height_array[1] = 16/2;
        }
    if (input_resolution <= INPUT_SIZE_720p_RANGE)
        me_context_ptr->hme_decimation = pcs_ptr->enc_mode <= ENC_MR ? ONE_DECIMATION_HME : TWO_DECIMATION_HME;
    else
        me_context_ptr->hme_decimation = TWO_DECIMATION_HME;

    // Scale up the MIN ME area if low frame rate
    uint8_t  low_frame_rate_flag = (scs_ptr->static_config.frame_rate >> 16) < 50 ? 1 : 0;
    if (low_frame_rate_flag) {
        me_context_ptr->search_area_width = (me_context_ptr->search_area_width * 3) / 2;
        me_context_ptr->search_area_height = (me_context_ptr->search_area_height * 3) / 2;
    }

    me_context_ptr->update_hme_search_center_flag = 1;

    if (input_resolution <= INPUT_SIZE_480p_RANGE)
        me_context_ptr->update_hme_search_center_flag = 0;

    return NULL;
};
void set_me_hme_ref_prune_ctrls(MeContext* context_ptr, uint8_t prune_level) {
    MeHmeRefPruneCtrls* me_hme_prune_ctrls = &context_ptr->me_hme_prune_ctrls;

    switch (prune_level)
    {
    case 0:
        me_hme_prune_ctrls->enable_me_hme_ref_pruning = 0;
        me_hme_prune_ctrls->prune_ref_if_hme_sad_dev_bigger_than_th = (uint16_t)~0;
        me_hme_prune_ctrls->prune_ref_if_me_sad_dev_bigger_than_th = (uint16_t)~0;
        break;
    case 1:
        me_hme_prune_ctrls->enable_me_hme_ref_pruning = 1;
        me_hme_prune_ctrls->prune_ref_if_hme_sad_dev_bigger_than_th = 160;
        me_hme_prune_ctrls->prune_ref_if_me_sad_dev_bigger_than_th = (uint16_t)~0;
        break;
    case 2:
        me_hme_prune_ctrls->enable_me_hme_ref_pruning = 1;
        me_hme_prune_ctrls->prune_ref_if_hme_sad_dev_bigger_than_th = 80;
        me_hme_prune_ctrls->prune_ref_if_me_sad_dev_bigger_than_th = 60;
        break;
    case 3:
        me_hme_prune_ctrls->enable_me_hme_ref_pruning = 1;
        me_hme_prune_ctrls->prune_ref_if_hme_sad_dev_bigger_than_th = 50;
        me_hme_prune_ctrls->prune_ref_if_me_sad_dev_bigger_than_th = 60;
        break;
    case 4:
        me_hme_prune_ctrls->enable_me_hme_ref_pruning = 1;
        me_hme_prune_ctrls->prune_ref_if_hme_sad_dev_bigger_than_th = 30;
        me_hme_prune_ctrls->prune_ref_if_me_sad_dev_bigger_than_th = 60;
        break;
    default:
        assert(0);
        break;
    }
}

void set_me_sr_adjustment_ctrls(MeContext* context_ptr, uint8_t sr_adjustment_level) {
    MeSrCtrls* me_sr_adjustment_ctrls = &context_ptr->me_sr_adjustment_ctrls;

    switch (sr_adjustment_level)
    {
    case 0:
        me_sr_adjustment_ctrls->enable_me_sr_adjustment = 0;
        break;
    case 1:
        me_sr_adjustment_ctrls->enable_me_sr_adjustment = 1;
        me_sr_adjustment_ctrls->reduce_me_sr_based_on_mv_length_th = 0;
        me_sr_adjustment_ctrls->stationary_hme_sad_abs_th = 100;
        me_sr_adjustment_ctrls->stationary_me_sr_divisor = 16;
        me_sr_adjustment_ctrls->reduce_me_sr_based_on_hme_sad_abs_th = 100;
        me_sr_adjustment_ctrls->me_sr_divisor_for_low_hme_sad = 8;
        break;
    case 2:
        me_sr_adjustment_ctrls->enable_me_sr_adjustment = 1;
        me_sr_adjustment_ctrls->reduce_me_sr_based_on_mv_length_th = 4;
        me_sr_adjustment_ctrls->stationary_hme_sad_abs_th = 12000;
        me_sr_adjustment_ctrls->stationary_me_sr_divisor = 8;
        me_sr_adjustment_ctrls->reduce_me_sr_based_on_hme_sad_abs_th = 6000;
        me_sr_adjustment_ctrls->me_sr_divisor_for_low_hme_sad = 8;
        break;
    default:
        assert(0);
        break;
    }
}
/******************************************************
* GM controls
******************************************************/
void set_gm_controls(PictureParentControlSet *pcs_ptr, uint8_t gm_level)
{
    GmControls *gm_ctrls = &pcs_ptr->gm_ctrls;
    switch (gm_level)
    {
    case 0:
        gm_ctrls->enabled = 0;
        break;
    case 1:
        gm_ctrls->enabled = 1;
        gm_ctrls->identiy_exit = 0;
        gm_ctrls->rotzoom_model_only = 0;
        gm_ctrls->bipred_only = 0;
        break;
    case 2:
        gm_ctrls->enabled = 1;
        gm_ctrls->identiy_exit = 1;
        gm_ctrls->rotzoom_model_only = 0;
        gm_ctrls->bipred_only = 0;
        break;
    case 3:
        gm_ctrls->enabled = 1;
        gm_ctrls->identiy_exit = 1;
        gm_ctrls->rotzoom_model_only = 1;
        gm_ctrls->bipred_only = 0;
        break;
    case 4:
        gm_ctrls->enabled = 1;
        gm_ctrls->identiy_exit = 1;
        gm_ctrls->rotzoom_model_only = 1;
        gm_ctrls->bipred_only = 1;
        break;
    default:
        assert(0);
        break;
    }
}
/******************************************************
* Derive ME Settings for OQ
  Input   : encoder mode and tune
  Output  : ME Kernel signal(s)
******************************************************/
EbErrorType signal_derivation_me_kernel_oq(SequenceControlSet *       scs_ptr,
                                           PictureParentControlSet *  pcs_ptr,
                                           MotionEstimationContext_t *context_ptr) {
    EbErrorType return_error = EB_ErrorNone;

    EbEncMode enc_mode = pcs_ptr->enc_mode;
    EbInputResolution input_resolution = scs_ptr->input_resolution;
    // Set ME/HME search regions
    if (scs_ptr->static_config.use_default_me_hme)
        set_me_hme_params_oq(
            context_ptr->me_context_ptr, pcs_ptr, scs_ptr, input_resolution);

    else
        set_me_hme_params_from_config(scs_ptr, context_ptr->me_context_ptr);

    // Set HME flags
    context_ptr->me_context_ptr->enable_hme_flag        = pcs_ptr->enable_hme_flag;
    context_ptr->me_context_ptr->enable_hme_level0_flag = pcs_ptr->enable_hme_level0_flag;
    context_ptr->me_context_ptr->enable_hme_level1_flag = pcs_ptr->enable_hme_level1_flag;
    context_ptr->me_context_ptr->enable_hme_level2_flag = pcs_ptr->enable_hme_level2_flag;
    // HME Search Method
    if (enc_mode <= ENC_MRS)
        context_ptr->me_context_ptr->hme_search_method = FULL_SAD_SEARCH;
    else
        context_ptr->me_context_ptr->hme_search_method = SUB_SAD_SEARCH;
    if (enc_mode <= ENC_MRS)
        context_ptr->me_context_ptr->me_search_method = FULL_SAD_SEARCH;
    else
        context_ptr->me_context_ptr->me_search_method = SUB_SAD_SEARCH;
    uint8_t gm_level = 0;
    if (scs_ptr->static_config.enable_global_motion == EB_TRUE &&
        pcs_ptr->frame_superres_enabled == EB_FALSE) {
        if (enc_mode <= ENC_M1)
            gm_level = 2;
        else if (enc_mode <= ENC_M6)
            gm_level = 3;
        else
            gm_level = pcs_ptr->is_used_as_reference_flag ? 4 : 0;
    }
    set_gm_controls(pcs_ptr, gm_level);

    // Set hme/me based reference pruning level (0-4)
    if (enc_mode <= ENC_MR)
            set_me_hme_ref_prune_ctrls(context_ptr->me_context_ptr, 0);
    else if (enc_mode <= ENC_M2)
            set_me_hme_ref_prune_ctrls(context_ptr->me_context_ptr, 2);
    else
            set_me_hme_ref_prune_ctrls(context_ptr->me_context_ptr, 4);
    // Set hme-based me sr adjustment level
    if (enc_mode <= ENC_MRS)
        set_me_sr_adjustment_ctrls(context_ptr->me_context_ptr, 0);
    else
        set_me_sr_adjustment_ctrls(context_ptr->me_context_ptr, 2);
    return return_error;
};
void open_loop_first_pass(struct PictureParentControlSet *ppcs_ptr,
                                 MotionEstimationContext_t *me_context_ptr, int32_t segment_index);

/******************************************************
* Derive ME Settings for first pass
  Input   : encoder mode and tune
  Output  : ME Kernel signal(s)
******************************************************/
EbErrorType first_pass_signal_derivation_me_kernel(
    SequenceControlSet        *scs_ptr,
    PictureParentControlSet   *pcs_ptr,
    MotionEstimationContext_t   *context_ptr) ;

/************************************************
 * Set ME/HME Params for Altref Temporal Filtering
 ************************************************/
void *tf_set_me_hme_params_oq(MeContext *me_context_ptr, PictureParentControlSet *pcs_ptr,
                              SequenceControlSet *scs_ptr, EbInputResolution input_resolution) {
    UNUSED(scs_ptr);
    UNUSED(pcs_ptr);
    // HME/ME default settings
    me_context_ptr->number_hme_search_region_in_width  = 2;
    me_context_ptr->number_hme_search_region_in_height = 2;
    // Set the minimum ME search area
    if (pcs_ptr->enc_mode <= ENC_M4) {
        me_context_ptr->search_area_width = me_context_ptr->search_area_height = input_resolution <= INPUT_SIZE_480p_RANGE ? 60 : 16;
        me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = input_resolution <= INPUT_SIZE_480p_RANGE ? 120 : 32;
}
    else {
        me_context_ptr->search_area_width = me_context_ptr->search_area_height = 8;
        me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 16;
    }
    if (pcs_ptr->enc_mode <= ENC_M4)
        me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = 30;
    else
        me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = 16;
    if (pcs_ptr->enc_mode <= ENC_M4)
        me_context_ptr->hme_level0_max_total_search_area_width = me_context_ptr->hme_level0_max_total_search_area_height = 60;
    else
        me_context_ptr->hme_level0_max_total_search_area_width = me_context_ptr->hme_level0_max_total_search_area_height = 32;


    me_context_ptr->hme_level0_max_search_area_in_width_array[0] =
        me_context_ptr->hme_level0_max_search_area_in_width_array[1] =
        me_context_ptr->hme_level0_max_total_search_area_width / me_context_ptr->number_hme_search_region_in_width;
    me_context_ptr->hme_level0_max_search_area_in_height_array[0] =
        me_context_ptr->hme_level0_max_search_area_in_height_array[1] =
        me_context_ptr->hme_level0_max_total_search_area_height / me_context_ptr->number_hme_search_region_in_height;
    me_context_ptr->hme_level0_search_area_in_width_array[0] =
        me_context_ptr->hme_level0_search_area_in_width_array[1] =
        me_context_ptr->hme_level0_total_search_area_width / me_context_ptr->number_hme_search_region_in_width;
    me_context_ptr->hme_level0_search_area_in_height_array[0] =
        me_context_ptr->hme_level0_search_area_in_height_array[1] =
        me_context_ptr->hme_level0_total_search_area_height / me_context_ptr->number_hme_search_region_in_height;

    me_context_ptr->hme_level1_search_area_in_width_array[0] =
        me_context_ptr->hme_level1_search_area_in_width_array[1] =
        me_context_ptr->hme_level1_search_area_in_height_array[0] =
        me_context_ptr->hme_level1_search_area_in_height_array[1] = 16;

    me_context_ptr->hme_level2_search_area_in_width_array[0] =
        me_context_ptr->hme_level2_search_area_in_width_array[1] =
        me_context_ptr->hme_level2_search_area_in_height_array[0] =
        me_context_ptr->hme_level2_search_area_in_height_array[1] = 16;

    me_context_ptr->hme_decimation = TWO_DECIMATION_HME;
    me_context_ptr->update_hme_search_center_flag = 1;

    if (input_resolution <= INPUT_SIZE_480p_RANGE)
        me_context_ptr->update_hme_search_center_flag = 0;
    return NULL;
};
/******************************************************
* Derive ME Settings for OQ for Altref Temporal Filtering
  Input   : encoder mode and tune
  Output  : ME Kernel signal(s)
******************************************************/
EbErrorType tf_signal_derivation_me_kernel_oq(SequenceControlSet *       scs_ptr,
                                              PictureParentControlSet *  pcs_ptr,
                                              MotionEstimationContext_t *context_ptr) {
    EbErrorType return_error = EB_ErrorNone;
    EbInputResolution input_resolution = scs_ptr->input_resolution;
    // Set ME/HME search regions
    tf_set_me_hme_params_oq(
        context_ptr->me_context_ptr, pcs_ptr, scs_ptr, input_resolution);
    // Set HME flags
    context_ptr->me_context_ptr->enable_hme_flag        = pcs_ptr->tf_enable_hme_flag;
    context_ptr->me_context_ptr->enable_hme_level0_flag = pcs_ptr->tf_enable_hme_level0_flag;
    context_ptr->me_context_ptr->enable_hme_level1_flag = pcs_ptr->tf_enable_hme_level1_flag;
    context_ptr->me_context_ptr->enable_hme_level2_flag = pcs_ptr->tf_enable_hme_level2_flag;
    // HME Search Method
        context_ptr->me_context_ptr->hme_search_method = FULL_SAD_SEARCH;
    // ME Search Method
        context_ptr->me_context_ptr->me_search_method = SUB_SAD_SEARCH;

    // Set hme/me based reference pruning level (0-4)
    // Ref pruning is disallowed for TF in motion_estimate_sb()
    set_me_hme_ref_prune_ctrls(context_ptr->me_context_ptr, 0);

    // Set hme-based me sr adjustment level
    // ME SR adjustment is disallowed for TF in motion_estimate_sb()
    set_me_sr_adjustment_ctrls(context_ptr->me_context_ptr, 0);
    return return_error;
};

static void motion_estimation_context_dctor(EbPtr p) {
    EbThreadContext *          thread_context_ptr = (EbThreadContext *)p;
    MotionEstimationContext_t *obj = (MotionEstimationContext_t *)thread_context_ptr->priv;
    EB_DELETE(obj->me_context_ptr);
    EB_FREE_ARRAY(obj);
}

/************************************************
 * Motion Analysis Context Constructor
 ************************************************/
EbErrorType motion_estimation_context_ctor(EbThreadContext *  thread_context_ptr,
                                           const EbEncHandle *enc_handle_ptr, int index) {
    MotionEstimationContext_t *context_ptr;

    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv  = context_ptr;
    thread_context_ptr->dctor = motion_estimation_context_dctor;
    context_ptr->picture_decision_results_input_fifo_ptr = svt_system_resource_get_consumer_fifo(
        enc_handle_ptr->picture_decision_results_resource_ptr, index);
    context_ptr->motion_estimation_results_output_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->motion_estimation_results_resource_ptr, index);
    EB_NEW(context_ptr->me_context_ptr,
        me_context_ctor);
    return EB_ErrorNone;
}

/***************************************************************************************************
* ZZ Decimated SAD Computation
***************************************************************************************************/
EbErrorType compute_decimated_zz_sad(MotionEstimationContext_t *context_ptr, PictureParentControlSet *pcs_ptr,
                                     EbPictureBufferDesc *sixteenth_decimated_picture_ptr,
                                     uint32_t x_sb_start_index, uint32_t x_sb_end_index,
                                     uint32_t y_sb_start_index, uint32_t y_sb_end_index) {
    EbErrorType return_error = EB_ErrorNone;

    PictureParentControlSet *previous_picture_control_set_wrapper_ptr =
        ((PictureParentControlSet *)pcs_ptr->previous_picture_control_set_wrapper_ptr->object_ptr);
    EbPictureBufferDesc *prev_input_picture_full =
        previous_picture_control_set_wrapper_ptr->enhanced_picture_ptr;

    uint32_t sb_index;

    uint32_t sb_width;
    uint32_t sb_height;

    uint32_t decimated_sb_width;
    uint32_t decimated_sb_height;

    uint32_t sb_origin_x;
    uint32_t sb_origin_y;

    uint32_t blk_displacement_decimated;
    uint32_t blk_displacement_full;

    uint32_t decimated_sb_collocated_sad;

    uint32_t x_sb_index;
    uint32_t y_sb_index;

    for (y_sb_index = y_sb_start_index; y_sb_index < y_sb_end_index; ++y_sb_index) {
        for (x_sb_index = x_sb_start_index; x_sb_index < x_sb_end_index; ++x_sb_index) {
            sb_index            = x_sb_index + y_sb_index * pcs_ptr->picture_sb_width;
            SbParams *sb_params = &pcs_ptr->sb_params_array[sb_index];

            sb_origin_x = sb_params->origin_x;
            sb_origin_y = sb_params->origin_y;

            sb_width  = sb_params->width;
            sb_height = sb_params->height;

            decimated_sb_width  = sb_width >> 2;
            decimated_sb_height = sb_height >> 2;

            decimated_sb_collocated_sad = 0;

            if (sb_params->is_complete_sb) {
                blk_displacement_decimated =
                    (sixteenth_decimated_picture_ptr->origin_y + (sb_origin_y >> 2)) *
                        sixteenth_decimated_picture_ptr->stride_y +
                    sixteenth_decimated_picture_ptr->origin_x + (sb_origin_x >> 2);
                blk_displacement_full = (prev_input_picture_full->origin_y + sb_origin_y) *
                                            prev_input_picture_full->stride_y +
                                        (prev_input_picture_full->origin_x + sb_origin_x);

                // 1/16 collocated SB decimation
                decimation_2d(&prev_input_picture_full->buffer_y[blk_displacement_full],
                              prev_input_picture_full->stride_y,
                              BLOCK_SIZE_64,
                              BLOCK_SIZE_64,
                              context_ptr->me_context_ptr->sixteenth_sb_buffer,
                              context_ptr->me_context_ptr->sixteenth_sb_buffer_stride,
                              4);

                // ZZ SAD between 1/16 current & 1/16 collocated
                decimated_sb_collocated_sad = svt_nxm_sad_kernel(
                    &(sixteenth_decimated_picture_ptr->buffer_y[blk_displacement_decimated]),
                    sixteenth_decimated_picture_ptr->stride_y,
                    context_ptr->me_context_ptr->sixteenth_sb_buffer,
                    context_ptr->me_context_ptr->sixteenth_sb_buffer_stride,
                    16,
                    16);
            } else {
                decimated_sb_collocated_sad = (uint32_t)~0;
            }
            // Keep track of non moving SBs for QP modulation
            if (decimated_sb_collocated_sad < ((decimated_sb_width * decimated_sb_height) * 2))
                previous_picture_control_set_wrapper_ptr->non_moving_index_array[sb_index] =
                    BEA_CLASS_0_ZZ_COST;
            else if (decimated_sb_collocated_sad < ((decimated_sb_width * decimated_sb_height) * 4))
                previous_picture_control_set_wrapper_ptr->non_moving_index_array[sb_index] =
                    BEA_CLASS_1_ZZ_COST;
            else if (decimated_sb_collocated_sad < ((decimated_sb_width * decimated_sb_height) * 8))
                previous_picture_control_set_wrapper_ptr->non_moving_index_array[sb_index] =
                    BEA_CLASS_2_ZZ_COST;
            else
                previous_picture_control_set_wrapper_ptr->non_moving_index_array[sb_index] =
                    BEA_CLASS_3_ZZ_COST;
        }
    }

    return return_error;
}
/***************************************************************************************************
* ZZ Decimated SSD Computation
***************************************************************************************************/
EbErrorType compute_zz_ssd(/*MotionEstimationContext_t *context_ptr, */PictureParentControlSet *pcs_ptr,
    uint32_t x_sb_start_index, uint32_t x_sb_end_index,
    uint32_t y_sb_start_index, uint32_t y_sb_end_index) {
    EbErrorType return_error = EB_ErrorNone;

    PictureParentControlSet *previous_picture_control_set_wrapper_ptr =
        ((PictureParentControlSet *)pcs_ptr->previous_picture_control_set_wrapper_ptr->object_ptr);
    EbPictureBufferDesc *prev_input_picture_full =
        previous_picture_control_set_wrapper_ptr->enhanced_picture_ptr;
    EbPictureBufferDesc *input_picture_ptr = pcs_ptr->enhanced_unscaled_picture_ptr;
    SequenceControlSet *     scs_ptr = pcs_ptr->scs_ptr;
    const uint32_t mb_cols = (scs_ptr->seq_header.max_frame_width + FORCED_BLK_SIZE - 1) / FORCED_BLK_SIZE;

    uint32_t sb_width;
    uint32_t sb_height;
    uint32_t blk_width;
    uint32_t blk_height;
    uint32_t sb_origin_x;
    uint32_t sb_origin_y;
    uint32_t blk_origin_x;
    uint32_t blk_origin_y;

    uint32_t input_origin_index;

    uint32_t x_sb_index;
    uint32_t y_sb_index;
    EbSpatialFullDistType spatial_full_dist_type_fun = svt_spatial_full_distortion_kernel;

    for (y_sb_index = y_sb_start_index; y_sb_index < y_sb_end_index; ++y_sb_index) {
        for (x_sb_index = x_sb_start_index; x_sb_index < x_sb_end_index; ++x_sb_index) {
            sb_origin_x = x_sb_index * scs_ptr->sb_sz;
            sb_origin_y = y_sb_index * scs_ptr->sb_sz;

            sb_width =
                (pcs_ptr->aligned_width - sb_origin_x) < BLOCK_SIZE_64
                ? pcs_ptr->aligned_width - sb_origin_x
                : BLOCK_SIZE_64;
            sb_height =
                (pcs_ptr->aligned_height - sb_origin_y) < BLOCK_SIZE_64
                ? pcs_ptr->aligned_height - sb_origin_y
                : BLOCK_SIZE_64;

            for (uint32_t y_blk_index = 0; y_blk_index < (sb_height + FORCED_BLK_SIZE - 1) / FORCED_BLK_SIZE; ++y_blk_index) {
                for (uint32_t x_blk_index = 0; x_blk_index < (sb_width + FORCED_BLK_SIZE-1) / FORCED_BLK_SIZE; ++x_blk_index) {
                    blk_origin_x = sb_origin_x + x_blk_index * FORCED_BLK_SIZE;
                    blk_origin_y = sb_origin_y + y_blk_index* FORCED_BLK_SIZE;

                    blk_width =
                        (pcs_ptr->aligned_width - blk_origin_x) < FORCED_BLK_SIZE
                        ? pcs_ptr->aligned_width - blk_origin_x
                        : FORCED_BLK_SIZE;
                    blk_height =
                        (pcs_ptr->aligned_height - blk_origin_y) < FORCED_BLK_SIZE
                        ? pcs_ptr->aligned_height - blk_origin_y
                        : FORCED_BLK_SIZE;

                    input_origin_index = (input_picture_ptr->origin_y + blk_origin_y) *
                        input_picture_ptr->stride_y +
                        (input_picture_ptr->origin_x + blk_origin_x);

                    pcs_ptr->firstpass_data.raw_motion_err_list[(blk_origin_y >> 4) * mb_cols + (blk_origin_x >> 4)]
                        = (uint32_t)( spatial_full_dist_type_fun(input_picture_ptr->buffer_y,
                            input_origin_index,
                            input_picture_ptr->stride_y,
                            prev_input_picture_full->buffer_y,
                            input_origin_index,
                            input_picture_ptr->stride_y,
                            blk_width,
                            blk_height));
                }
            }
        }
    }

    return return_error;
}

/************************************************
 * Motion Analysis Kernel
 * The Motion Analysis performs  Motion Estimation
 * This process has access to the current input picture as well as
 * the input pictures, which the current picture references according
 * to the prediction structure pattern.  The Motion Analysis process is multithreaded,
 * so pictures can be processed out of order as long as all inputs are available.
 ************************************************/
void *motion_estimation_kernel(void *input_ptr) {
    EbThreadContext *          thread_context_ptr = (EbThreadContext *)input_ptr;
    MotionEstimationContext_t *context_ptr = (MotionEstimationContext_t *)thread_context_ptr->priv;
    EbObjectWrapper *          in_results_wrapper_ptr;
    EbObjectWrapper *          out_results_wrapper_ptr;
    for (;;) {
        // Get Input Full Object
        EB_GET_FULL_OBJECT(context_ptr->picture_decision_results_input_fifo_ptr,
                           &in_results_wrapper_ptr);
        PictureDecisionResults *in_results_ptr = (PictureDecisionResults *)
                                                     in_results_wrapper_ptr->object_ptr;
        PictureParentControlSet *pcs_ptr = (PictureParentControlSet *)
                                               in_results_ptr->pcs_wrapper_ptr->object_ptr;
        SequenceControlSet * scs_ptr = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;
        context_ptr->me_context_ptr->me_type =
            in_results_ptr->task_type == 1 ? ME_MCTF :
            in_results_ptr->task_type == 0 ? ME_OPEN_LOOP : ME_FIRST_PASS;

        // Lambda Assignement
        if (scs_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS) {
            if (pcs_ptr->temporal_layer_index == 0)
                context_ptr->me_context_ptr->lambda =
                    lambda_mode_decision_ra_sad[pcs_ptr->picture_qp];
            else if (pcs_ptr->temporal_layer_index < 3)
                context_ptr->me_context_ptr->lambda =
                    lambda_mode_decision_ra_sad_qp_scaling_l1[pcs_ptr->picture_qp];
            else
                context_ptr->me_context_ptr->lambda =
                    lambda_mode_decision_ra_sad_qp_scaling_l3[pcs_ptr->picture_qp];
        } else {
            if (pcs_ptr->temporal_layer_index == 0)
                context_ptr->me_context_ptr->lambda =
                    lambda_mode_decision_ld_sad[pcs_ptr->picture_qp];
            else
                context_ptr->me_context_ptr->lambda =
                    lambda_mode_decision_ld_sad_qp_scaling[pcs_ptr->picture_qp];
        }
        if (in_results_ptr->task_type == 0) {
            // ME Kernel Signal(s) derivation
            if (use_output_stat(scs_ptr))
                first_pass_signal_derivation_me_kernel(scs_ptr, pcs_ptr, context_ptr);
            else
                signal_derivation_me_kernel_oq(scs_ptr, pcs_ptr, context_ptr);
            EbPictureBufferDesc *sixteenth_picture_ptr = NULL;
            EbPictureBufferDesc *quarter_picture_ptr = NULL;
            EbPictureBufferDesc *input_padded_picture_ptr = NULL;
            EbPictureBufferDesc *input_picture_ptr = NULL;
            EbPaReferenceObject *pa_ref_obj_ = NULL;
            if (!scs_ptr->in_loop_me) {
                pa_ref_obj_ = (EbPaReferenceObject *)pcs_ptr->pa_reference_picture_wrapper_ptr->object_ptr;
                // Set 1/4 and 1/16 ME input buffer(s); filtered or decimated
                quarter_picture_ptr =
                    (scs_ptr->down_sampling_method_me_search == ME_FILTERED_DOWNSAMPLED)
                    ? (EbPictureBufferDesc *)pa_ref_obj_->quarter_filtered_picture_ptr
                    : (EbPictureBufferDesc *)pa_ref_obj_->quarter_decimated_picture_ptr;

                sixteenth_picture_ptr =
                    (scs_ptr->down_sampling_method_me_search == ME_FILTERED_DOWNSAMPLED)
                    ? (EbPictureBufferDesc *)pa_ref_obj_->sixteenth_filtered_picture_ptr
                    : (EbPictureBufferDesc *)pa_ref_obj_->sixteenth_decimated_picture_ptr;
                input_padded_picture_ptr = (EbPictureBufferDesc *)pa_ref_obj_->input_padded_picture_ptr;
            }
            input_picture_ptr = pcs_ptr->enhanced_unscaled_picture_ptr;
            // Segments
            uint32_t segment_index   = in_results_ptr->segment_index;
            uint32_t pic_width_in_sb = (pcs_ptr->aligned_width + scs_ptr->sb_sz - 1) /
                scs_ptr->sb_sz;
            uint32_t picture_height_in_sb = (pcs_ptr->aligned_height + scs_ptr->sb_sz - 1) /
                scs_ptr->sb_sz;
            uint32_t y_segment_index;
            uint32_t x_segment_index;
            SEGMENT_CONVERT_IDX_TO_XY(
                segment_index, x_segment_index, y_segment_index, pcs_ptr->me_segments_column_count);
            uint32_t x_sb_start_index = SEGMENT_START_IDX(
                x_segment_index, pic_width_in_sb, pcs_ptr->me_segments_column_count);
            uint32_t x_sb_end_index = SEGMENT_END_IDX(
                x_segment_index, pic_width_in_sb, pcs_ptr->me_segments_column_count);
            uint32_t y_sb_start_index = SEGMENT_START_IDX(
                y_segment_index, picture_height_in_sb, pcs_ptr->me_segments_row_count);
            uint32_t y_sb_end_index = SEGMENT_END_IDX(
                y_segment_index, picture_height_in_sb, pcs_ptr->me_segments_row_count);
            EbBool skip_me = EB_FALSE;
            if (use_output_stat(scs_ptr))
                skip_me = EB_TRUE;
            // skip me for the first pass. ME is already performed
            if (!skip_me) {
            // *** MOTION ESTIMATION CODE ***
            if (pcs_ptr->slice_type != I_SLICE && !scs_ptr->in_loop_me) {
                // Use scaled source references if resolution of the reference is different that of the input
                use_scaled_source_refs_if_needed(pcs_ptr,
                                                 input_picture_ptr,
                                                 pa_ref_obj_,
                                                 &input_padded_picture_ptr,
                                                 &quarter_picture_ptr,
                                                 &sixteenth_picture_ptr);

                // SB Loop
                for (uint32_t y_sb_index = y_sb_start_index; y_sb_index < y_sb_end_index;
                     ++y_sb_index) {
                    for (uint32_t x_sb_index = x_sb_start_index; x_sb_index < x_sb_end_index;
                         ++x_sb_index) {
                        uint32_t sb_index = (uint16_t)(x_sb_index + y_sb_index * pic_width_in_sb);
                        uint32_t sb_origin_x = x_sb_index * scs_ptr->sb_sz;
                        uint32_t sb_origin_y = y_sb_index * scs_ptr->sb_sz;

                        uint32_t sb_width = (pcs_ptr->aligned_width - sb_origin_x) < BLOCK_SIZE_64
                            ? pcs_ptr->aligned_width - sb_origin_x
                            : BLOCK_SIZE_64;

                        // Load the SB from the input to the intermediate SB buffer
                        uint32_t buffer_index = (input_picture_ptr->origin_y + sb_origin_y) *
                                input_picture_ptr->stride_y +
                            input_picture_ptr->origin_x + sb_origin_x;
                        for (unsigned sb_row = 0; sb_row < BLOCK_SIZE_64; sb_row++) {
                            svt_memcpy(
                                &(context_ptr->me_context_ptr->sb_buffer[sb_row * BLOCK_SIZE_64]),
                                &(input_picture_ptr
                                      ->buffer_y[buffer_index +
                                                 sb_row * input_picture_ptr->stride_y]),
                                sizeof(uint8_t) * BLOCK_SIZE_64);
                        }
#ifdef ARCH_X86_64
                        uint8_t *src_ptr   = &input_padded_picture_ptr->buffer_y[buffer_index];
                        uint32_t sb_height = (pcs_ptr->aligned_height - sb_origin_y) < BLOCK_SIZE_64
                            ? pcs_ptr->aligned_height - sb_origin_y
                            : BLOCK_SIZE_64;
                        //_MM_HINT_T0     //_MM_HINT_T1    //_MM_HINT_T2//_MM_HINT_NTA
                        for (uint32_t i = 0; i < sb_height; i++) {
                            char const *p = (char const *)(src_ptr +
                                                           i * input_padded_picture_ptr->stride_y);

                            _mm_prefetch(p, _MM_HINT_T2);
                        }
#endif

                        context_ptr->me_context_ptr->sb_src_ptr =
                            &input_padded_picture_ptr->buffer_y[buffer_index];
                        context_ptr->me_context_ptr->sb_src_stride =
                            input_padded_picture_ptr->stride_y;
                        // Load the 1/4 decimated SB from the 1/4 decimated input to the 1/4 intermediate SB buffer
                        if (context_ptr->me_context_ptr->enable_hme_level1_flag) {
                            buffer_index = (quarter_picture_ptr->origin_y + (sb_origin_y >> 1)) *
                                    quarter_picture_ptr->stride_y +
                                quarter_picture_ptr->origin_x + (sb_origin_x >> 1);

                            for (unsigned sb_row = 0; sb_row < (BLOCK_SIZE_64 >> 1); sb_row++) {
                                svt_memcpy(
                                    &(context_ptr->me_context_ptr->quarter_sb_buffer
                                          [sb_row *
                                           context_ptr->me_context_ptr->quarter_sb_buffer_stride]),
                                    &(quarter_picture_ptr
                                          ->buffer_y[buffer_index +
                                                     sb_row * quarter_picture_ptr->stride_y]),
                                    sizeof(uint8_t) * (sb_width >> 1));
                            }
                        }

                        // Load the 1/16 decimated SB from the 1/16 decimated input to the 1/16 intermediate SB buffer
                        if (context_ptr->me_context_ptr->enable_hme_level0_flag) {
                            buffer_index = (sixteenth_picture_ptr->origin_y + (sb_origin_y >> 2)) *
                                    sixteenth_picture_ptr->stride_y +
                                sixteenth_picture_ptr->origin_x + (sb_origin_x >> 2);

                            uint8_t *frame_ptr = &sixteenth_picture_ptr->buffer_y[buffer_index];
                            uint8_t *local_ptr = context_ptr->me_context_ptr->sixteenth_sb_buffer;
                            for (unsigned sb_row = 0; sb_row < (BLOCK_SIZE_64 >> 2);
                                 sb_row += context_ptr->me_context_ptr->hme_search_method ==
                                         FULL_SAD_SEARCH
                                     ? 1
                                     : 2) {
                                svt_memcpy(local_ptr, frame_ptr, (sb_width >> 2) * sizeof(uint8_t));
                                local_ptr += 16;
                                frame_ptr += sixteenth_picture_ptr->stride_y
                                    << (context_ptr->me_context_ptr->hme_search_method !=
                                        FULL_SAD_SEARCH);
                            }
                        }
                        context_ptr->me_context_ptr->me_type = ME_OPEN_LOOP;
                        context_ptr->me_context_ptr->num_of_list_to_search =
                            (pcs_ptr->slice_type == P_SLICE) ? (uint32_t)REF_LIST_0 : (uint32_t)REF_LIST_1;

                        context_ptr->me_context_ptr->num_of_ref_pic_to_search[0] = pcs_ptr->ref_list0_count_try;
                        if (pcs_ptr->slice_type == B_SLICE)
                            context_ptr->me_context_ptr->num_of_ref_pic_to_search[1] = pcs_ptr->ref_list1_count_try;
                        context_ptr->me_context_ptr->temporal_layer_index = pcs_ptr->temporal_layer_index;
                        context_ptr->me_context_ptr->is_used_as_reference_flag = pcs_ptr->is_used_as_reference_flag;

                        for (int i = 0; i<= context_ptr->me_context_ptr->num_of_list_to_search; i++) {
                            for (int j=0; j< context_ptr->me_context_ptr->num_of_ref_pic_to_search[i];j++) {
                                EbPaReferenceObject* reference_object =
                                    (EbPaReferenceObject *)pcs_ptr->ref_pa_pic_ptr_array[i][j]->object_ptr;
                                context_ptr->me_context_ptr->me_ds_ref_array[i][j].picture_ptr =
                                    reference_object->input_padded_picture_ptr;
                                if (scs_ptr->down_sampling_method_me_search == ME_FILTERED_DOWNSAMPLED) {
                                    context_ptr->me_context_ptr->me_ds_ref_array[i][j].quarter_picture_ptr =
                                        reference_object->quarter_filtered_picture_ptr;
                                    context_ptr->me_context_ptr->me_ds_ref_array[i][j].sixteenth_picture_ptr =
                                        reference_object->sixteenth_filtered_picture_ptr;
                                } else {
                                    context_ptr->me_context_ptr->me_ds_ref_array[i][j].quarter_picture_ptr =
                                        reference_object->quarter_decimated_picture_ptr;
                                    context_ptr->me_context_ptr->me_ds_ref_array[i][j].sixteenth_picture_ptr =
                                        reference_object->sixteenth_decimated_picture_ptr;
                                }
                                context_ptr->me_context_ptr->me_ds_ref_array[i][j].picture_number = reference_object->picture_number;
                            }
                        }

                        motion_estimate_sb(pcs_ptr,
                                           sb_index,
                                           sb_origin_x,
                                           sb_origin_y,
                                           context_ptr->me_context_ptr,
                                           input_picture_ptr);
                        svt_block_on_mutex(pcs_ptr->me_processed_sb_mutex);
                        pcs_ptr->me_processed_sb_count++;
                        // We need to finish ME for all SBs to do GM
                        if (pcs_ptr->me_processed_sb_count == pcs_ptr->sb_total_count) {
                            if (pcs_ptr->gm_ctrls.enabled)
                                global_motion_estimation(

                                    pcs_ptr, input_picture_ptr);
                            else
                            // Initilize global motion to be OFF when GM is OFF
                                memset(pcs_ptr->is_global_motion, EB_FALSE, MAX_NUM_OF_REF_PIC_LIST * REF_LIST_MAX_DEPTH);
                        }

                        svt_release_mutex(pcs_ptr->me_processed_sb_mutex);

                    }
                }
            }
            if (
                scs_ptr->in_loop_ois == 0 &&
                (!scs_ptr->in_loop_me || pcs_ptr->slice_type == I_SLICE) &&
                scs_ptr->static_config.enable_tpl_la)
                for (uint32_t y_sb_index = y_sb_start_index; y_sb_index < y_sb_end_index;
                     ++y_sb_index)
                    for (uint32_t x_sb_index = x_sb_start_index; x_sb_index < x_sb_end_index;
                         ++x_sb_index) {
                        uint32_t sb_index = (uint16_t)(x_sb_index + y_sb_index * pic_width_in_sb);
                        open_loop_intra_search_mb(pcs_ptr, sb_index, input_picture_ptr);
                    }
            // ZZ SADs Computation
            // 1 lookahead frame is needed to get valid (0,0) SAD
            if (scs_ptr->static_config.look_ahead_distance != 0 &&
                    pcs_ptr->picture_number > 0 &&
                    !scs_ptr->in_loop_me)
                // when DG is ON, the ZZ SADs are computed @ the PD process
                // ZZ SADs Computation using decimated picture
                compute_decimated_zz_sad(
                    context_ptr,
                    pcs_ptr,
                    (EbPictureBufferDesc *)pa_ref_obj_
                        ->sixteenth_decimated_picture_ptr, // Hsan: always use decimated for ZZ SAD derivation until studying the trade offs and regenerating the activity threshold
                    x_sb_start_index,
                    x_sb_end_index,
                    y_sb_start_index,
                    y_sb_end_index);

            if (scs_ptr->static_config.look_ahead_distance != 0 &&
                    pcs_ptr->picture_number > 0 &&
                    scs_ptr->in_loop_me)
                compute_decimated_zz_sad(
                    context_ptr,
                    pcs_ptr,
                    pcs_ptr->ds_pics.sixteenth_picture_ptr,
                    x_sb_start_index,
                    x_sb_end_index,
                    y_sb_start_index,
                    y_sb_end_index);

            if (scs_ptr->static_config.rate_control_mode && !use_input_stat(scs_ptr) && !scs_ptr->lap_enabled) {
                // Calculate the ME Distortion and OIS Historgrams
                svt_block_on_mutex(pcs_ptr->rc_distortion_histogram_mutex);

                if (scs_ptr->static_config.rate_control_mode
                    && !(use_input_stat(scs_ptr) && scs_ptr->static_config.rate_control_mode == 1) //skip 2pass VBR
                    ) {
                    for (uint32_t y_sb_index = y_sb_start_index; y_sb_index < y_sb_end_index;
                        ++y_sb_index)
                        for (uint32_t x_sb_index = x_sb_start_index; x_sb_index < x_sb_end_index;
                            ++x_sb_index) {
                        uint32_t sb_origin_x = x_sb_index * scs_ptr->sb_sz;
                        uint32_t sb_origin_y = y_sb_index * scs_ptr->sb_sz;
                        uint32_t sb_width = (pcs_ptr->aligned_width - sb_origin_x) <
                            BLOCK_SIZE_64
                            ? pcs_ptr->aligned_width - sb_origin_x
                            : BLOCK_SIZE_64;
                        uint32_t sb_height = (pcs_ptr->aligned_height - sb_origin_y) <
                            BLOCK_SIZE_64
                            ? pcs_ptr->aligned_height - sb_origin_y
                            : BLOCK_SIZE_64;

                        uint32_t sb_index = (uint16_t)(x_sb_index +
                            y_sb_index * pic_width_in_sb);
                        pcs_ptr->inter_sad_interval_index[sb_index] = 0;
                        pcs_ptr->intra_sad_interval_index[sb_index] = 0;

                        if (sb_width == BLOCK_SIZE_64 && sb_height == BLOCK_SIZE_64) {
                            if (pcs_ptr->slice_type != I_SLICE && !scs_ptr->in_loop_me) {
                                uint16_t sad_interval_index = (uint16_t)(
                                    pcs_ptr->rc_me_distortion[sb_index] >>
                                    (12 - SAD_PRECISION_INTERVAL)); //change 12 to 2*log2(64)

                                sad_interval_index = (uint16_t)(sad_interval_index >> 2);
                                if (sad_interval_index > (NUMBER_OF_SAD_INTERVALS >> 1) - 1) {
                                    uint16_t sad_interval_index_temp = sad_interval_index -
                                        ((NUMBER_OF_SAD_INTERVALS >> 1) - 1);

                                    sad_interval_index = ((NUMBER_OF_SAD_INTERVALS >> 1) - 1) +
                                        (sad_interval_index_temp >> 3);
                                }
                                if (sad_interval_index >= NUMBER_OF_SAD_INTERVALS - 1)
                                    sad_interval_index = NUMBER_OF_SAD_INTERVALS - 1;

                                pcs_ptr->inter_sad_interval_index[sb_index] = sad_interval_index;

                                pcs_ptr->me_distortion_histogram[sad_interval_index]++;
                            }

                            uint32_t intra_sad_interval_index =
                                pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_64x64] >> 4;
                            intra_sad_interval_index = (uint16_t)(intra_sad_interval_index >>
                                2);
                            if (intra_sad_interval_index > (NUMBER_OF_SAD_INTERVALS >> 1) - 1) {
                                uint32_t sad_interval_index_temp = intra_sad_interval_index -
                                    ((NUMBER_OF_SAD_INTERVALS >> 1) - 1);

                                intra_sad_interval_index = ((NUMBER_OF_SAD_INTERVALS >> 1) -
                                    1) +
                                    (sad_interval_index_temp >> 3);
                            }
                            if (intra_sad_interval_index >= NUMBER_OF_SAD_INTERVALS - 1)
                                intra_sad_interval_index = NUMBER_OF_SAD_INTERVALS - 1;

                            pcs_ptr->intra_sad_interval_index[sb_index] =
                                intra_sad_interval_index;

                            pcs_ptr->ois_distortion_histogram[intra_sad_interval_index]++;

                            ++pcs_ptr->full_sb_count;
                        }
                    }
                }

                svt_release_mutex(pcs_ptr->rc_distortion_histogram_mutex);
            }
            }
            // Get Empty Results Object
            svt_get_empty_object(context_ptr->motion_estimation_results_output_fifo_ptr,
                                &out_results_wrapper_ptr);

            MotionEstimationResults *out_results_ptr = (MotionEstimationResults *)
                                                           out_results_wrapper_ptr->object_ptr;
            out_results_ptr->pcs_wrapper_ptr = in_results_ptr->pcs_wrapper_ptr;
            out_results_ptr->segment_index   = segment_index;

            // Release the Input Results
            svt_release_object(in_results_wrapper_ptr);

            // Post the Full Results Object
            svt_post_full_object(out_results_wrapper_ptr);
        } else if (in_results_ptr->task_type == 1) {
            // ME Kernel Signal(s) derivation
            tf_signal_derivation_me_kernel_oq(scs_ptr, pcs_ptr, context_ptr);

            // temporal filtering start
            context_ptr->me_context_ptr->me_type = ME_MCTF;
            svt_av1_init_temporal_filtering(
                pcs_ptr->temp_filt_pcs_list, pcs_ptr, context_ptr, in_results_ptr->segment_index);

            // Release the Input Results
            svt_release_object(in_results_wrapper_ptr);
        }
        else {
            // ME Kernel Signal(s) derivation
            first_pass_signal_derivation_me_kernel(scs_ptr, pcs_ptr, context_ptr);

            //For first pass compute_decimated_zz_sad() is skipped, and non_moving_index_array[] become uninitialized
            init_zz_cost_info((PictureParentControlSet *)
                                  pcs_ptr->previous_picture_control_set_wrapper_ptr->object_ptr);

            // first pass start
            context_ptr->me_context_ptr->me_type = ME_FIRST_PASS;
            open_loop_first_pass(
                pcs_ptr, context_ptr, in_results_ptr->segment_index);

            // Release the Input Results
            svt_release_object(in_results_wrapper_ptr);
        }
    }

    return NULL;
}
// inloop ME ctor
EbErrorType ime_context_ctor(EbThreadContext *  thread_context_ptr,
        const EbEncHandle *enc_handle_ptr, int index) {
    InLoopMeContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv = context_ptr;
    thread_context_ptr->dctor = motion_estimation_context_dctor;
    context_ptr->input_fifo_ptr = svt_system_resource_get_consumer_fifo(
            enc_handle_ptr->pic_mgr_res_srm, index);

    context_ptr->output_fifo_ptr = svt_system_resource_get_producer_fifo(
            enc_handle_ptr->rate_control_tasks_resource_ptr, index);

    EB_NEW(context_ptr->me_context_ptr,
            me_context_ctor);
    return EB_ErrorNone;
}
// Lambda Assignement
static void init_lambda(InLoopMeContext *context_ptr,
    SequenceControlSet *scs_ptr,
    PictureParentControlSet *ppcs_ptr) {
    uint8_t temporal_layer_index = context_ptr->me_context_ptr->temporal_layer_index;
    if (scs_ptr->static_config.pred_structure == EB_PRED_RANDOM_ACCESS) {
        if (temporal_layer_index == 0)
            context_ptr->me_context_ptr->lambda =
                lambda_mode_decision_ra_sad[ppcs_ptr->picture_qp];
        else if (temporal_layer_index < 3)
            context_ptr->me_context_ptr->lambda =
                lambda_mode_decision_ra_sad_qp_scaling_l1[ppcs_ptr->picture_qp];
        else
            context_ptr->me_context_ptr->lambda =
                lambda_mode_decision_ra_sad_qp_scaling_l3[ppcs_ptr->picture_qp];
    } else {
        if (temporal_layer_index == 0)
            context_ptr->me_context_ptr->lambda =
                lambda_mode_decision_ld_sad[ppcs_ptr->picture_qp];
        else
            context_ptr->me_context_ptr->lambda =
                lambda_mode_decision_ld_sad_qp_scaling[ppcs_ptr->picture_qp];
    }
}

// Get ME buffer sb based
static void prepare_sb_me_buffer(InLoopMeContext *context_ptr,
    PictureParentControlSet *ppcs_ptr,
    uint32_t sb_origin_x, uint32_t sb_origin_y) {

    // Get 1/4 and 1/16 ME reference buffer(s); filtered or decimated
    EbDownScaledObject *src_ds_object =
        (EbDownScaledObject*)
        ppcs_ptr->down_scaled_picture_wrapper_ptr->object_ptr;
    EbPictureBufferDesc *quarter_picture_ptr = src_ds_object->quarter_picture_ptr;
    EbPictureBufferDesc *sixteenth_picture_ptr = src_ds_object->sixteenth_picture_ptr;
    EbPictureBufferDesc *input_picture_ptr = ppcs_ptr->enhanced_picture_ptr;

    uint32_t sb_width =
        (ppcs_ptr->aligned_width - sb_origin_x) < BLOCK_SIZE_64
        ? ppcs_ptr->aligned_width - sb_origin_x
        : BLOCK_SIZE_64;

    // Load the SB from the input to the intermediate SB buffer
    uint32_t buffer_index = (input_picture_ptr->origin_y + sb_origin_y) *
        input_picture_ptr->stride_y +
        input_picture_ptr->origin_x + sb_origin_x;

#ifdef ARCH_X86
    uint8_t *src_ptr = &input_picture_ptr->buffer_y[buffer_index];
    uint32_t sb_height =
        (ppcs_ptr->aligned_height - sb_origin_y) < BLOCK_SIZE_64
        ? ppcs_ptr->aligned_height  - sb_origin_y
        : BLOCK_SIZE_64;

    //_MM_HINT_T0     //_MM_HINT_T1    //_MM_HINT_T2//_MM_HINT_NTA
    uint32_t i;
    for (i = 0; i < sb_height; i++) {
        char const *p =
            (char const *)(src_ptr +
                    i * input_picture_ptr->stride_y);
        _mm_prefetch(p, _MM_HINT_T2);
    }
#endif
    context_ptr->me_context_ptr->sb_src_ptr =
        &input_picture_ptr->buffer_y[buffer_index];
    context_ptr->me_context_ptr->sb_src_stride =
        input_picture_ptr->stride_y;

    // Load the 1/4 decimated SB from the 1/4 decimated input to the 1/4 intermediate SB buffer
    if (context_ptr->me_context_ptr->enable_hme_level1_flag) {
        buffer_index = (quarter_picture_ptr->origin_y + (sb_origin_y >> 1)) *
            quarter_picture_ptr->stride_y +
            quarter_picture_ptr->origin_x + (sb_origin_x >> 1);
        for (uint32_t sb_row = 0; sb_row < (BLOCK_SIZE_64 >> 1); sb_row++) {
            EB_MEMCPY(
                    (&(context_ptr->me_context_ptr
                       ->quarter_sb_buffer[sb_row *
                       context_ptr->me_context_ptr
                       ->quarter_sb_buffer_stride])),
                    (&(quarter_picture_ptr
                       ->buffer_y[buffer_index +
                       sb_row * quarter_picture_ptr->stride_y])),
                    (sb_width >> 1) * sizeof(uint8_t));
        }
    }

    // Load the 1/16 decimated SB from the 1/16 decimated input to the 1/16 intermediate SB buffer
    if (context_ptr->me_context_ptr->enable_hme_level0_flag) {
        buffer_index = (sixteenth_picture_ptr->origin_y + (sb_origin_y >> 2)) *
            sixteenth_picture_ptr->stride_y +
            sixteenth_picture_ptr->origin_x + (sb_origin_x >> 2);

        uint8_t *frame_ptr = &sixteenth_picture_ptr->buffer_y[buffer_index];
        uint8_t *local_ptr =
            context_ptr->me_context_ptr->sixteenth_sb_buffer;
        if (context_ptr->me_context_ptr->hme_search_method ==
                FULL_SAD_SEARCH) {
            for (uint32_t sb_row = 0; sb_row < (BLOCK_SIZE_64 >> 2); sb_row++) {
                EB_MEMCPY(local_ptr,
                        frame_ptr,
                        (sb_width >> 2) * sizeof(uint8_t));
                local_ptr += 16;
                frame_ptr += sixteenth_picture_ptr->stride_y;
            }
        } else {
            for (uint32_t sb_row = 0; sb_row < (BLOCK_SIZE_64 >> 2); sb_row += 2) {
                EB_MEMCPY(local_ptr,
                        frame_ptr,
                        (sb_width >> 2) * sizeof(uint8_t));
                local_ptr += 16;
                frame_ptr += sixteenth_picture_ptr->stride_y << 1;
            }
        }
    }
}

/************************************************
 * inloop Motion Analysis Kernel
 * The Motion Analysis performs  Motion Estimation
 * This process has access to the current input picture as well as
 * the reference pictures, which can be input or reconstructed
 ************************************************/
void *inloop_me_kernel(void *input_ptr) {
    EbThreadContext *          thread_context_ptr = (EbThreadContext *)input_ptr;
    InLoopMeContext *context_ptr = (InLoopMeContext *)thread_context_ptr->priv;


    EbObjectWrapper *       in_results_wrapper_ptr;
    PictureManagerResults   *in_results_ptr;

    EbObjectWrapper *        out_results_wrapper_ptr;

    EbPictureBufferDesc *input_picture_ptr;

    uint32_t pic_width_in_sb = 0;
    uint32_t picture_height_in_sb;
    uint32_t sb_origin_x;
    uint32_t sb_origin_y;

    // Segments
    uint32_t segment_index = 0;
    uint32_t x_segment_index;
    uint32_t y_segment_index;
    uint32_t x_sb_start_index = 0;
    uint32_t x_sb_end_index = 0;
    uint32_t y_sb_start_index = 0;
    uint32_t y_sb_end_index = 0;
    uint32_t segment_col_count = 0;
    uint32_t segment_row_count = 0;

    EbBool skip_me = EB_FALSE;
    for (;;) {
        // Get Input Full Object
        EB_GET_FULL_OBJECT(context_ptr->input_fifo_ptr,
                &in_results_wrapper_ptr);

        in_results_ptr = (PictureManagerResults *)in_results_wrapper_ptr->object_ptr;
        PictureParentControlSet* ppcs_ptr = (PictureParentControlSet*)in_results_ptr->pcs_wrapper_ptr->object_ptr;
        SequenceControlSet* scs_ptr =
            (SequenceControlSet *)ppcs_ptr->scs_wrapper_ptr->object_ptr;
        uint8_t task_type = in_results_ptr->task_type;

        // iME get ppcs input, and output pcs to RC kernel
        if (scs_ptr->in_loop_me) {
            input_picture_ptr = ppcs_ptr->enhanced_picture_ptr;

            segment_col_count = ppcs_ptr->inloop_me_segments_column_count;
            segment_row_count = ppcs_ptr->inloop_me_segments_row_count;

            if (task_type != 0) {
                signal_derivation_me_kernel_oq(scs_ptr, ppcs_ptr, (MotionEstimationContext_t*)context_ptr);

                // TPL ME
                segment_col_count = ppcs_ptr->tpl_me_segments_column_count;
                segment_row_count = ppcs_ptr->tpl_me_segments_row_count;
                context_ptr->me_context_ptr->me_type = ME_TPL;
                context_ptr->me_context_ptr->num_of_list_to_search = (in_results_ptr->tpl_ref_list1_count > 0) ?
                    REF_LIST_1 : REF_LIST_0;
                context_ptr->me_context_ptr->num_of_ref_pic_to_search[0] = in_results_ptr->tpl_ref_list0_count;
                context_ptr->me_context_ptr->num_of_ref_pic_to_search[1] = in_results_ptr->tpl_ref_list1_count;
                context_ptr->me_context_ptr->temporal_layer_index = in_results_ptr->temporal_layer_index;
                context_ptr->me_context_ptr->is_used_as_reference_flag = in_results_ptr->is_used_as_reference_flag;
                for (int i = 0; i<= context_ptr->me_context_ptr->num_of_list_to_search; i++) {
                    for (int j=0; j< context_ptr->me_context_ptr->num_of_ref_pic_to_search[i];j++) {
                        context_ptr->me_context_ptr->me_ds_ref_array[i][j] =
                            ppcs_ptr->tpl_data.tpl_ref_ds_ptr_array[i][j];
                    }
                }
                skip_me = EB_FALSE;
            } else if (ppcs_ptr->slice_type != I_SLICE) {
                // ME Kernel Signal(s) derivation
                signal_derivation_me_kernel_oq(scs_ptr, ppcs_ptr, (MotionEstimationContext_t*)context_ptr);

                context_ptr->me_context_ptr->me_type = ME_CLOSE_LOOP;
                context_ptr->me_context_ptr->num_of_list_to_search =
                    (ppcs_ptr->slice_type == P_SLICE) ? (uint32_t)REF_LIST_0 : (uint32_t)REF_LIST_1;
                context_ptr->me_context_ptr->num_of_ref_pic_to_search[0] = ppcs_ptr->ref_list0_count_try;
                if (ppcs_ptr->slice_type == B_SLICE)
                    context_ptr->me_context_ptr->num_of_ref_pic_to_search[1] = ppcs_ptr->ref_list1_count_try;
                context_ptr->me_context_ptr->temporal_layer_index = ppcs_ptr->temporal_layer_index;
                context_ptr->me_context_ptr->is_used_as_reference_flag = ppcs_ptr->is_used_as_reference_flag;

                for (int i = 0; i<= context_ptr->me_context_ptr->num_of_list_to_search; i++) {
                    for (int j=0; j< context_ptr->me_context_ptr->num_of_ref_pic_to_search[i];j++) {
                        EbReferenceObject* inl_reference_object =
                            ppcs_ptr->child_pcs->ref_pic_ptr_array[i][j]->object_ptr;

                        context_ptr->me_context_ptr->me_ds_ref_array[i][j].picture_ptr =
                            inl_reference_object->reference_picture;
                        context_ptr->me_context_ptr->me_ds_ref_array[i][j].sixteenth_picture_ptr =
                            inl_reference_object->sixteenth_reference_picture;
                        context_ptr->me_context_ptr->me_ds_ref_array[i][j].quarter_picture_ptr =
                            inl_reference_object->quarter_reference_picture;
                        context_ptr->me_context_ptr->me_ds_ref_array[i][j].picture_number =
                            ppcs_ptr->ref_pic_poc_array[i][j];

                        context_ptr->me_context_ptr->me_ds_ref_array[i][j].picture_ptr =
                            inl_reference_object->input_picture;
                        context_ptr->me_context_ptr->me_ds_ref_array[i][j].sixteenth_picture_ptr =
                            inl_reference_object->sixteenth_input_picture;
                        context_ptr->me_context_ptr->me_ds_ref_array[i][j].quarter_picture_ptr =
                            inl_reference_object->quarter_input_picture;
                    }
                }
                skip_me = ppcs_ptr->tpl_me_done;
            }
            if(use_output_stat(scs_ptr))
                skip_me = EB_TRUE;
            // Segments
            segment_index = in_results_ptr->segment_index;

            // case TPL ON , trailing frames ON
            // 1. trailing pic, task_type == 2 (tplME) ,skip_me == 0 , slice_type not set
            // 2. non-trailing pic, task_type == 2 (tplME) ,skip_me == 0, slice_type set
            // 3. non-trailing pic, task_type == 0 (iME) ,skip_me == 1, slice_type set
            if (!skip_me && (ppcs_ptr->slice_type != I_SLICE || task_type != 0)) {
                // Lambda Assignement
                init_lambda(context_ptr, scs_ptr, ppcs_ptr);

                pic_width_in_sb =
                    (ppcs_ptr->aligned_width + scs_ptr->sb_sz - 1) / scs_ptr->sb_sz;
                picture_height_in_sb =
                    (ppcs_ptr->aligned_height + scs_ptr->sb_sz - 1) / scs_ptr->sb_sz;
                SEGMENT_CONVERT_IDX_TO_XY(
                        segment_index, x_segment_index, y_segment_index, segment_col_count);
                x_sb_start_index = SEGMENT_START_IDX(
                        x_segment_index, pic_width_in_sb, segment_col_count);
                x_sb_end_index = SEGMENT_END_IDX(
                        x_segment_index, pic_width_in_sb, segment_col_count);
                y_sb_start_index = SEGMENT_START_IDX(
                        y_segment_index, picture_height_in_sb, segment_row_count);
                y_sb_end_index = SEGMENT_END_IDX(
                        y_segment_index, picture_height_in_sb, segment_row_count);

                // SB Loop
                for (uint32_t y_sb_index = y_sb_start_index; y_sb_index < y_sb_end_index; ++y_sb_index) {
                    for (uint32_t x_sb_index = x_sb_start_index; x_sb_index < x_sb_end_index; ++x_sb_index) {
                        uint32_t sb_index = x_sb_index + y_sb_index * pic_width_in_sb;
                        sb_origin_x = x_sb_index * scs_ptr->sb_sz;
                        sb_origin_y = y_sb_index * scs_ptr->sb_sz;
                        prepare_sb_me_buffer(context_ptr, ppcs_ptr, sb_origin_x, sb_origin_y);

                        motion_estimate_sb(ppcs_ptr,
                                sb_index,
                                sb_origin_x,
                                sb_origin_y,
                                context_ptr->me_context_ptr,
                                input_picture_ptr);
                        {
                            svt_block_on_mutex(ppcs_ptr->me_processed_sb_mutex);
                            ppcs_ptr->me_processed_sb_count++;
                            svt_release_mutex(ppcs_ptr->me_processed_sb_mutex);
                        }
                    }
                }
            }
            if (task_type == 0) {


                if (scs_ptr->static_config.enable_tpl_la) {
                    if (segment_index == 0) {
                        if (ppcs_ptr->gm_ctrls.enabled && ppcs_ptr->slice_type != I_SLICE)
                            global_motion_estimation_inl(
                                ppcs_ptr, input_picture_ptr);
                        else
                            // Initilize global motion to be OFF for all references frames.
                            memset(ppcs_ptr->is_global_motion, EB_FALSE, MAX_NUM_OF_REF_PIC_LIST * REF_LIST_MAX_DEPTH);
                    }

                }
                else {
                    svt_block_on_mutex(ppcs_ptr->me_processed_sb_mutex);
                    if (ppcs_ptr->me_processed_sb_count == ppcs_ptr->sb_total_count) {
                        if (ppcs_ptr->gm_ctrls.enabled && ppcs_ptr->slice_type != I_SLICE)
                            global_motion_estimation_inl(
                                ppcs_ptr, input_picture_ptr);
                        else
                            // Initilize global motion to be OFF for all references frames.
                            memset(ppcs_ptr->is_global_motion, EB_FALSE, MAX_NUM_OF_REF_PIC_LIST * REF_LIST_MAX_DEPTH);
                    }
                    svt_release_mutex(ppcs_ptr->me_processed_sb_mutex);

                }

                svt_get_empty_object(context_ptr->output_fifo_ptr,
                        &out_results_wrapper_ptr);

                RateControlTasks * rate_control_tasks_ptr = (RateControlTasks *)out_results_wrapper_ptr->object_ptr;
                rate_control_tasks_ptr->pcs_wrapper_ptr = ppcs_ptr->child_pcs->c_pcs_wrapper_ptr;
                rate_control_tasks_ptr->task_type = RC_INPUT;
                rate_control_tasks_ptr->segment_index = segment_index;

                // Release the Input Results
                svt_release_object(in_results_wrapper_ptr);

                // Post the Full Results Object
                svt_post_full_object(out_results_wrapper_ptr);
            } else {
                // TPL ME
                // Doing OIS search for TPL
                if (scs_ptr->in_loop_ois == 0)
                if (scs_ptr->static_config.enable_tpl_la) {
                    for (uint32_t y_sb_index = y_sb_start_index; y_sb_index < y_sb_end_index;
                            ++y_sb_index) {
                        for (uint32_t x_sb_index = x_sb_start_index; x_sb_index < x_sb_end_index;
                                ++x_sb_index) {
                            uint32_t sb_index = x_sb_index + y_sb_index * pic_width_in_sb;
                            open_loop_intra_search_mb(ppcs_ptr, sb_index, input_picture_ptr);
                        }
                    }
                }
                svt_block_on_mutex(ppcs_ptr->tpl_me_mutex);
                ppcs_ptr->tpl_me_seg_acc++;

                if (ppcs_ptr->tpl_me_seg_acc ==
                        ppcs_ptr->tpl_me_segments_total_count) {

                    svt_post_semaphore(ppcs_ptr->tpl_me_done_semaphore);
                }
                svt_release_mutex(ppcs_ptr->tpl_me_mutex);
                svt_release_object(in_results_wrapper_ptr);
            }
        } else {
            // Dummy in loop kernel
            svt_get_empty_object(context_ptr->output_fifo_ptr,
                    &out_results_wrapper_ptr);

            RateControlTasks * rate_control_tasks_ptr = (RateControlTasks *)out_results_wrapper_ptr->object_ptr;
            rate_control_tasks_ptr->pcs_wrapper_ptr = ppcs_ptr->child_pcs->c_pcs_wrapper_ptr;
            rate_control_tasks_ptr->task_type = RC_INPUT;
            rate_control_tasks_ptr->segment_index = segment_index;

            // Release the Input Results
            svt_release_object(in_results_wrapper_ptr);

            // Post the Full Results Object
            svt_post_full_object(out_results_wrapper_ptr);
        }
    }

    return NULL;
}
