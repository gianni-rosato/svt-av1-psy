/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
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
#ifdef ARCH_X86
#include "emmintrin.h"
#endif
#include "EbTemporalFiltering.h"
#include "EbGlobalMotionEstimation.h"

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
    uint8_t hme_me_level =
        scs_ptr->use_output_stat_file ? pcs_ptr->snd_pass_enc_mode : pcs_ptr->enc_mode;
#if !CS2_ADOPTIONS_1
    if (hme_me_level <= ENC_M1) hme_me_level = ENC_M0;
#endif
    // HME/ME default settings
    me_context_ptr->number_hme_search_region_in_width  = 2;
    me_context_ptr->number_hme_search_region_in_height = 2;

    uint8_t sc_content_detected = pcs_ptr->sc_content_detected;

#if DIST_BASED_ME_SEARCH_AREA
    // sequence frame rate based mutiplier
    uint8_t  fr_rate_mult_num   = sc_content_detected ? 1 :
        (scs_ptr->static_config.frame_rate >> 16) < 50 ? 3 : 1;
    uint8_t  fr_rate_mult_denum = sc_content_detected ? 1 :
        (scs_ptr->static_config.frame_rate >> 16) < 50 ? 2 : 1;

#endif

    // HME Level0
    me_context_ptr->hme_level0_total_search_area_width =
        hme_level0_total_search_area_width[sc_content_detected][input_resolution][hme_me_level];
    me_context_ptr->hme_level0_total_search_area_height =
        hme_level0_total_search_area_height[sc_content_detected][input_resolution][hme_me_level];
    me_context_ptr->hme_level0_search_area_in_width_array[0] =
        hme_level0_search_area_in_width_array_right[sc_content_detected][input_resolution]
                                                   [hme_me_level];
    me_context_ptr->hme_level0_search_area_in_width_array[1] =
        hme_level0_search_area_in_width_array_left[sc_content_detected][input_resolution]
                                                  [hme_me_level];
    me_context_ptr->hme_level0_search_area_in_height_array[0] =
        hme_level0_search_area_in_height_array_top[sc_content_detected][input_resolution]
                                                  [hme_me_level];
    me_context_ptr->hme_level0_search_area_in_height_array[1] =
        hme_level0_search_area_in_height_array_bottom[sc_content_detected][input_resolution]
                                                     [hme_me_level];
    // HME Level1
    me_context_ptr->hme_level1_search_area_in_width_array[0] =
        hme_level1_search_area_in_width_array_right[sc_content_detected][input_resolution]
                                                   [hme_me_level];
    me_context_ptr->hme_level1_search_area_in_width_array[1] =
        hme_level1_search_area_in_width_array_left[sc_content_detected][input_resolution]
                                                  [hme_me_level];
    me_context_ptr->hme_level1_search_area_in_height_array[0] =
        hme_level1_search_area_in_height_array_top[sc_content_detected][input_resolution]
                                                  [hme_me_level];
    me_context_ptr->hme_level1_search_area_in_height_array[1] =
        hme_level1_search_area_in_height_array_bottom[sc_content_detected][input_resolution]
                                                     [hme_me_level];
    // HME Level2
    me_context_ptr->hme_level2_search_area_in_width_array[0] =
        hme_level2_search_area_in_width_array_right[sc_content_detected][input_resolution]
                                                   [hme_me_level];
    me_context_ptr->hme_level2_search_area_in_width_array[1] =
        hme_level2_search_area_in_width_array_left[sc_content_detected][input_resolution]
                                                  [hme_me_level];
    me_context_ptr->hme_level2_search_area_in_height_array[0] =
        hme_level2_search_area_in_height_array_top[sc_content_detected][input_resolution]
                                                  [hme_me_level];
    me_context_ptr->hme_level2_search_area_in_height_array[1] =
        hme_level2_search_area_in_height_array_bottom[sc_content_detected][input_resolution]
                                                     [hme_me_level];

#if DIST_BASED_ME_SEARCH_AREA
    me_context_ptr->search_area_width  =
        (min_me_search_width[sc_content_detected][input_resolution][hme_me_level]
            *fr_rate_mult_num)/fr_rate_mult_denum ;
    me_context_ptr->search_area_height =
        (min_me_search_height[sc_content_detected][input_resolution][hme_me_level]
            *fr_rate_mult_num)/fr_rate_mult_denum ;
#else
    // ME
    me_context_ptr->search_area_width =
        search_area_width[sc_content_detected][input_resolution][hme_me_level];
    me_context_ptr->search_area_height =
        search_area_height[sc_content_detected][input_resolution][hme_me_level];
#endif
    assert(me_context_ptr->search_area_width <= MAX_SEARCH_AREA_WIDTH &&
           "increase MAX_SEARCH_AREA_WIDTH");
    assert(me_context_ptr->search_area_height <= MAX_SEARCH_AREA_HEIGHT &&
           "increase MAX_SEARCH_AREA_HEIGHT");

    me_context_ptr->update_hme_search_center_flag = 1;

    if (input_resolution <= INPUT_SIZE_576p_RANGE_OR_LOWER)
        me_context_ptr->update_hme_search_center_flag = 0;

    return NULL;
};

/******************************************************
* Derive ME Settings for OQ
  Input   : encoder mode and tune
  Output  : ME Kernel signal(s)
******************************************************/
EbErrorType signal_derivation_me_kernel_oq(SequenceControlSet *       scs_ptr,
                                           PictureParentControlSet *  pcs_ptr,
                                           MotionEstimationContext_t *context_ptr) {
    EbErrorType return_error = EB_ErrorNone;

    uint8_t enc_mode =
        scs_ptr->use_output_stat_file ? pcs_ptr->snd_pass_enc_mode : pcs_ptr->enc_mode;
    EbInputResolution input_resolution = scs_ptr->input_resolution;
    uint8_t sc_content_detected = pcs_ptr->sc_content_detected;
#if DIST_BASED_ME_SEARCH_AREA
    uint8_t  hme_me_level = scs_ptr->use_output_stat_file ?
        pcs_ptr->snd_pass_enc_mode : pcs_ptr->enc_mode;

#if !CS2_ADOPTIONS_1
    if (hme_me_level <= ENC_M2)
        hme_me_level = ENC_M0;
#endif
#endif
    // Set ME/HME search regions
    if (scs_ptr->static_config.use_default_me_hme)
        set_me_hme_params_oq(
            context_ptr->me_context_ptr, pcs_ptr, scs_ptr, input_resolution);

    else
        set_me_hme_params_from_config(scs_ptr, context_ptr->me_context_ptr);

#if DIST_BASED_ME_SEARCH_AREA
    // ME
    context_ptr->me_context_ptr->max_me_search_width =
        max_me_search_width[sc_content_detected][input_resolution][hme_me_level];
    context_ptr->me_context_ptr->max_me_search_height =
        max_me_search_height[sc_content_detected][input_resolution][hme_me_level];
#endif
    if (sc_content_detected)
        context_ptr->me_context_ptr->fractional_search_method =
#if CS2_ADOPTIONS_1
            (enc_mode <= ENC_M1) ? FULL_SAD_SEARCH : SUB_SAD_SEARCH;
#else
            (enc_mode == ENC_M0) ? FULL_SAD_SEARCH : SUB_SAD_SEARCH;
#endif
    else if (enc_mode <= ENC_M6)
        context_ptr->me_context_ptr->fractional_search_method = SSD_SEARCH;
    else
        context_ptr->me_context_ptr->fractional_search_method = FULL_SAD_SEARCH;

    // Set HME flags
    context_ptr->me_context_ptr->enable_hme_flag        = pcs_ptr->enable_hme_flag;
    context_ptr->me_context_ptr->enable_hme_level0_flag = pcs_ptr->enable_hme_level0_flag;
    context_ptr->me_context_ptr->enable_hme_level1_flag = pcs_ptr->enable_hme_level1_flag;
    context_ptr->me_context_ptr->enable_hme_level2_flag = pcs_ptr->enable_hme_level2_flag;

    if (scs_ptr->static_config.enable_subpel == DEFAULT)
        // Set the default settings of subpel
        if (sc_content_detected)
#if !CS2_ADOPTIONS_1
            if (enc_mode <= ENC_M5)
                context_ptr->me_context_ptr->use_subpel_flag = 1;
            else
#endif
                context_ptr->me_context_ptr->use_subpel_flag = 0;
        else
            context_ptr->me_context_ptr->use_subpel_flag = 1;
    else
        context_ptr->me_context_ptr->use_subpel_flag = scs_ptr->static_config.enable_subpel;
#if CS2_ADOPTIONS_1
    if (enc_mode <= ENC_M0) {
        context_ptr->me_context_ptr->half_pel_mode =
            (sc_content_detected) ? REFINEMENT_HP_MODE : EX_HP_MODE;
#else
    if (MR_MODE) {
        context_ptr->me_context_ptr->half_pel_mode    = EX_HP_MODE;
    } else if (enc_mode <= ENC_M0) {
        context_ptr->me_context_ptr->half_pel_mode =
            (sc_content_detected) ? REFINEMENT_HP_MODE : EX_HP_MODE;
#endif
#if SWITCHED_HALF_PEL_MODE
    }else if (enc_mode <= ENC_M2) {
        context_ptr->me_context_ptr->half_pel_mode =
            (sc_content_detected) ? REFINEMENT_HP_MODE : SWITCHABLE_HP_MODE;
#endif
    } else {
        context_ptr->me_context_ptr->half_pel_mode    = REFINEMENT_HP_MODE;
    }

#if CS2_ADOPTIONS_1
    context_ptr->me_context_ptr->h_pel_search_wind = H_PEL_SEARCH_WIND_2;
#endif
    // Set fractional search model
    // 0: search all blocks
    // 1: selective based on Full-Search SAD & MV.
    // 2: off
    if (context_ptr->me_context_ptr->use_subpel_flag == 1) {
        if (enc_mode <= ENC_M6)
            context_ptr->me_context_ptr->fractional_search_model = 0;
        else
            context_ptr->me_context_ptr->fractional_search_model = 1;
    } else
        context_ptr->me_context_ptr->fractional_search_model = 2;
    // HME Search Method
#if CS2_ADOPTIONS_1
        context_ptr->me_context_ptr->hme_search_method = SUB_SAD_SEARCH;
#else
    if (sc_content_detected)
        if (enc_mode <= ENC_M6)
            context_ptr->me_context_ptr->hme_search_method = FULL_SAD_SEARCH;
        else
            context_ptr->me_context_ptr->hme_search_method = SUB_SAD_SEARCH;
    else
        context_ptr->me_context_ptr->hme_search_method = SUB_SAD_SEARCH;
#endif
    // ME Search Method
#if CS2_ADOPTIONS_1
        context_ptr->me_context_ptr->me_search_method = SUB_SAD_SEARCH;
#else
    if (sc_content_detected)
        if (enc_mode <= ENC_M5)
            context_ptr->me_context_ptr->me_search_method = FULL_SAD_SEARCH;
        else
            context_ptr->me_context_ptr->me_search_method = SUB_SAD_SEARCH;
    else
        context_ptr->me_context_ptr->me_search_method = SUB_SAD_SEARCH;
#endif

    if (scs_ptr->static_config.enable_global_motion == EB_TRUE) {
        if (enc_mode <= ENC_M1)
            context_ptr->me_context_ptr->compute_global_motion = EB_TRUE;
        else
            context_ptr->me_context_ptr->compute_global_motion = EB_FALSE;
    } else
        context_ptr->me_context_ptr->compute_global_motion = EB_FALSE;

    // Me nsq search levels.
    // 0: feature off -> perform nsq_search.
    // 1: perform me nsq_search only for the best refrenece picture.
    // 2: perform me nsq_search only for the nearest refrenece pictures.
    // 3: me nsq_search off.
    if (MR_MODE && pcs_ptr->sc_content_detected == 0)
        context_ptr->me_context_ptr->inherit_rec_mv_from_sq_block = 0;
    else
        context_ptr->me_context_ptr->inherit_rec_mv_from_sq_block = 2;

    return return_error;
};

/************************************************
 * Set ME/HME Params for Altref Temporal Filtering
 ************************************************/
void *tf_set_me_hme_params_oq(MeContext *me_context_ptr, PictureParentControlSet *pcs_ptr,
                              SequenceControlSet *scs_ptr, EbInputResolution input_resolution) {
    UNUSED(scs_ptr);
    uint8_t hme_me_level =
        scs_ptr->use_output_stat_file ? pcs_ptr->snd_pass_enc_mode : pcs_ptr->enc_mode;
    // HME/ME default settings
    me_context_ptr->number_hme_search_region_in_width  = 2;
    me_context_ptr->number_hme_search_region_in_height = 2;

    uint8_t sc_content_detected = pcs_ptr->sc_content_detected;

    // HME Level0
    me_context_ptr->hme_level0_total_search_area_width =
        tf_hme_level0_total_search_area_width[sc_content_detected][input_resolution][hme_me_level];
    me_context_ptr->hme_level0_total_search_area_height =
        tf_hme_level0_total_search_area_height[sc_content_detected][input_resolution][hme_me_level];
    me_context_ptr->hme_level0_search_area_in_width_array[0] =
        tf_hme_level0_search_area_in_width_array_right[sc_content_detected][input_resolution]
                                                      [hme_me_level];
    me_context_ptr->hme_level0_search_area_in_width_array[1] =
        tf_hme_level0_search_area_in_width_array_left[sc_content_detected][input_resolution]
                                                     [hme_me_level];
    me_context_ptr->hme_level0_search_area_in_height_array[0] =
        tf_hme_level0_search_area_in_height_array_top[sc_content_detected][input_resolution]
                                                     [hme_me_level];
    me_context_ptr->hme_level0_search_area_in_height_array[1] =
        tf_hme_level0_search_area_in_height_array_bottom[sc_content_detected][input_resolution]
                                                        [hme_me_level];
    // HME Level1
    me_context_ptr->hme_level1_search_area_in_width_array[0] =
        tf_hme_level1_search_area_in_width_array_right[sc_content_detected][input_resolution]
                                                      [hme_me_level];
    me_context_ptr->hme_level1_search_area_in_width_array[1] =
        tf_hme_level1_search_area_in_width_array_left[sc_content_detected][input_resolution]
                                                     [hme_me_level];
    me_context_ptr->hme_level1_search_area_in_height_array[0] =
        tf_hme_level1_search_area_in_height_array_top[sc_content_detected][input_resolution]
                                                     [hme_me_level];
    me_context_ptr->hme_level1_search_area_in_height_array[1] =
        tf_hme_level1_search_area_in_height_array_bottom[sc_content_detected][input_resolution]
                                                        [hme_me_level];
    // HME Level2
    me_context_ptr->hme_level2_search_area_in_width_array[0] =
        tf_hme_level2_search_area_in_width_array_right[sc_content_detected][input_resolution]
                                                      [hme_me_level];
    me_context_ptr->hme_level2_search_area_in_width_array[1] =
        tf_hme_level2_search_area_in_width_array_left[sc_content_detected][input_resolution]
                                                     [hme_me_level];
    me_context_ptr->hme_level2_search_area_in_height_array[0] =
        tf_hme_level2_search_area_in_height_array_top[sc_content_detected][input_resolution]
                                                     [hme_me_level];
    me_context_ptr->hme_level2_search_area_in_height_array[1] =
        tf_hme_level2_search_area_in_height_array_bottom[sc_content_detected][input_resolution]
                                                        [hme_me_level];

    // ME
#if DIST_BASED_ME_SEARCH_AREA
    me_context_ptr->search_area_width  =
        min_metf_search_width[sc_content_detected][input_resolution][hme_me_level];
    me_context_ptr->search_area_height =
        min_metf_search_height[sc_content_detected][input_resolution][hme_me_level];
#else
    me_context_ptr->search_area_width =
        tf_search_area_width[sc_content_detected][input_resolution][hme_me_level];
    me_context_ptr->search_area_height =
        tf_search_area_height[sc_content_detected][input_resolution][hme_me_level];
#endif
    assert(me_context_ptr->search_area_width <= MAX_SEARCH_AREA_WIDTH &&
           "increase MAX_SEARCH_AREA_WIDTH");
    assert(me_context_ptr->search_area_height <= MAX_SEARCH_AREA_HEIGHT &&
           "increase MAX_SEARCH_AREA_HEIGHT");

    me_context_ptr->update_hme_search_center_flag = 1;

    if (input_resolution <= INPUT_SIZE_576p_RANGE_OR_LOWER)
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
    uint8_t     enc_mode =
        scs_ptr->use_output_stat_file ? pcs_ptr->snd_pass_enc_mode : pcs_ptr->enc_mode;
    EbInputResolution input_resolution = scs_ptr->input_resolution;
    uint8_t sc_content_detected = pcs_ptr->sc_content_detected;
#if DIST_BASED_ME_SEARCH_AREA
    uint8_t  hme_me_level = scs_ptr->use_output_stat_file ?
        pcs_ptr->snd_pass_enc_mode : pcs_ptr->enc_mode;
#if !CS2_ADOPTIONS_1
    if (hme_me_level <= ENC_M2)
        hme_me_level = ENC_M0;
#endif

#endif
    // Set ME/HME search regions
    tf_set_me_hme_params_oq(
        context_ptr->me_context_ptr, pcs_ptr, scs_ptr, input_resolution);

#if DIST_BASED_ME_SEARCH_AREA // TF
    context_ptr->me_context_ptr->max_me_search_width =
        max_metf_search_width[sc_content_detected][input_resolution][hme_me_level];
    context_ptr->me_context_ptr->max_me_search_height =
        max_metf_search_height[sc_content_detected][input_resolution][hme_me_level];
#endif
    if (sc_content_detected)
        if (enc_mode <= ENC_M1)
#if CS2_ADOPTIONS_1
            context_ptr->me_context_ptr->fractional_search_method = FULL_SAD_SEARCH;
#else
            context_ptr->me_context_ptr->fractional_search_method =
                (enc_mode == ENC_M0) ? FULL_SAD_SEARCH : SSD_SEARCH;
#endif
        else
            context_ptr->me_context_ptr->fractional_search_method = SUB_SAD_SEARCH;
    else if (enc_mode <= ENC_M6)
        context_ptr->me_context_ptr->fractional_search_method = SSD_SEARCH;
    else
        context_ptr->me_context_ptr->fractional_search_method = FULL_SAD_SEARCH;

    // Set HME flags
    context_ptr->me_context_ptr->enable_hme_flag        = pcs_ptr->tf_enable_hme_flag;
    context_ptr->me_context_ptr->enable_hme_level0_flag = pcs_ptr->tf_enable_hme_level0_flag;
    context_ptr->me_context_ptr->enable_hme_level1_flag = pcs_ptr->tf_enable_hme_level1_flag;
    context_ptr->me_context_ptr->enable_hme_level2_flag = pcs_ptr->tf_enable_hme_level2_flag;
    if (scs_ptr->static_config.enable_subpel == DEFAULT)
        // Set the default settings of subpel
        if (sc_content_detected)
#if !CS2_ADOPTIONS_1
            if (enc_mode <= ENC_M1)
                context_ptr->me_context_ptr->use_subpel_flag = 1;
            else
#endif
                context_ptr->me_context_ptr->use_subpel_flag = 0;
        else
            context_ptr->me_context_ptr->use_subpel_flag = 1;
    else
        context_ptr->me_context_ptr->use_subpel_flag = scs_ptr->static_config.enable_subpel;
#if CS2_ADOPTIONS_1
    if (enc_mode <= ENC_M0) {
        context_ptr->me_context_ptr->half_pel_mode =
            (sc_content_detected) ? REFINEMENT_HP_MODE : EX_HP_MODE;
#else
    if (MR_MODE) {
        context_ptr->me_context_ptr->half_pel_mode    = EX_HP_MODE;
    } else if (enc_mode <= ENC_M0) {
        context_ptr->me_context_ptr->half_pel_mode =
            (sc_content_detected) ? REFINEMENT_HP_MODE : EX_HP_MODE;
#endif
#if SWITCHED_HALF_PEL_MODE
    }else if (enc_mode <= ENC_M1) {
        context_ptr->me_context_ptr->half_pel_mode =
            (sc_content_detected) ? REFINEMENT_HP_MODE : SWITCHABLE_HP_MODE;
#endif
    } else {
        context_ptr->me_context_ptr->half_pel_mode    = REFINEMENT_HP_MODE;
    }
#if CS2_ADOPTIONS_1
    context_ptr->me_context_ptr->h_pel_search_wind =   H_PEL_SEARCH_WIND_3;
#endif
    // Set fractional search model
    // 0: search all blocks
    // 1: selective based on Full-Search SAD & MV.
    // 2: off
    if (context_ptr->me_context_ptr->use_subpel_flag == 1) {
        if (enc_mode <= ENC_M6)
            context_ptr->me_context_ptr->fractional_search_model = 0;
        else
            context_ptr->me_context_ptr->fractional_search_model = 1;
    } else
        context_ptr->me_context_ptr->fractional_search_model = 2;

    // HME Search Method
    if (sc_content_detected)
        if (enc_mode <= ENC_M6)
            context_ptr->me_context_ptr->hme_search_method = FULL_SAD_SEARCH;
        else
            context_ptr->me_context_ptr->hme_search_method = SUB_SAD_SEARCH;
    else
        context_ptr->me_context_ptr->hme_search_method = FULL_SAD_SEARCH;
    // ME Search Method
    if (sc_content_detected)
        if (enc_mode <= ENC_M3)
            context_ptr->me_context_ptr->me_search_method = FULL_SAD_SEARCH;
        else
            context_ptr->me_context_ptr->me_search_method = SUB_SAD_SEARCH;
    else
        context_ptr->me_context_ptr->me_search_method =
            (enc_mode <= ENC_M4) ? FULL_SAD_SEARCH : SUB_SAD_SEARCH;
    // Me nsq search levels.
    // 0: feature off -> perform nsq_search.
    // 1: perform me nsq_search for the best refrenece picture.
    // 2: perform me nsq_search for the nearest refrenece pictures.
    // 3: me nsq_search off.
    context_ptr->me_context_ptr->inherit_rec_mv_from_sq_block = 0;

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
    const SequenceControlSet * scs_ptr;

    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv  = context_ptr;
    thread_context_ptr->dctor = motion_estimation_context_dctor;

    scs_ptr = enc_handle_ptr->scs_instance_array[0]->scs_ptr;

    context_ptr->picture_decision_results_input_fifo_ptr = eb_system_resource_get_consumer_fifo(
        enc_handle_ptr->picture_decision_results_resource_ptr, index);
    context_ptr->motion_estimation_results_output_fifo_ptr = eb_system_resource_get_producer_fifo(
        enc_handle_ptr->motion_estimation_results_resource_ptr, index);
    EB_NEW(context_ptr->me_context_ptr,
           me_context_ctor,
           scs_ptr->max_input_luma_width,
           scs_ptr->max_input_luma_height,
           scs_ptr->nsq_present,
           scs_ptr->mrp_mode);
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

            sb_width  = sb_params->width;
            sb_height = sb_params->height;

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
                decimated_sb_collocated_sad = nxm_sad_kernel(
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

    PictureParentControlSet *pcs_ptr;
    SequenceControlSet *     scs_ptr;

    EbObjectWrapper *       in_results_wrapper_ptr;
    PictureDecisionResults *in_results_ptr;

    EbObjectWrapper *        out_results_wrapper_ptr;
    MotionEstimationResults *out_results_ptr;

    EbPictureBufferDesc *input_picture_ptr;

    EbPictureBufferDesc *input_padded_picture_ptr;

    uint32_t buffer_index;

    uint32_t sb_index;
    uint32_t x_sb_index;
    uint32_t y_sb_index;
    uint32_t pic_width_in_sb;
    uint32_t picture_height_in_sb;
    uint32_t sb_origin_x;
    uint32_t sb_origin_y;
    uint32_t sb_width;
    uint32_t sb_height;
    uint32_t sb_row;

    EbPaReferenceObject *pa_ref_obj_;
    EbPictureBufferDesc *quarter_picture_ptr;
    EbPictureBufferDesc *sixteenth_picture_ptr;
    // Segments
    uint32_t segment_index;
    uint32_t x_segment_index;
    uint32_t y_segment_index;
    uint32_t x_sb_start_index;
    uint32_t x_sb_end_index;
    uint32_t y_sb_start_index;
    uint32_t y_sb_end_index;

    uint32_t intra_sad_interval_index;

    for (;;) {
        // Get Input Full Object
        EB_GET_FULL_OBJECT(context_ptr->picture_decision_results_input_fifo_ptr,
                           &in_results_wrapper_ptr);

        in_results_ptr = (PictureDecisionResults *)in_results_wrapper_ptr->object_ptr;
        pcs_ptr        = (PictureParentControlSet *)in_results_ptr->pcs_wrapper_ptr->object_ptr;
        scs_ptr        = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;

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

        input_picture_ptr = pcs_ptr->enhanced_unscaled_picture_ptr;

        context_ptr->me_context_ptr->me_alt_ref =
            in_results_ptr->task_type == 1 ? EB_TRUE : EB_FALSE;

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
            signal_derivation_me_kernel_oq(scs_ptr, pcs_ptr, context_ptr);

#if GLOBAL_WARPED_MOTION
            // Global motion estimation
            // Compute only for the first fragment.
            // TODO: create an other kernel ?
#if GLOBAL_WARPED_MOTION
            if (pcs_ptr->gm_level == GM_FULL || pcs_ptr->gm_level == GM_DOWN) {
#endif
                if (context_ptr->me_context_ptr->compute_global_motion &&
                    in_results_ptr->segment_index == 0)
                    global_motion_estimation(
                        pcs_ptr, context_ptr->me_context_ptr, input_picture_ptr);
#if GLOBAL_WARPED_MOTION
            }
#endif
#endif

            // Segments
            segment_index = in_results_ptr->segment_index;
            pic_width_in_sb =
                (pcs_ptr->aligned_width + scs_ptr->sb_sz - 1) / scs_ptr->sb_sz;
            picture_height_in_sb =
                (pcs_ptr->aligned_height + scs_ptr->sb_sz - 1) / scs_ptr->sb_sz;
            SEGMENT_CONVERT_IDX_TO_XY(
                segment_index, x_segment_index, y_segment_index, pcs_ptr->me_segments_column_count);
            x_sb_start_index = SEGMENT_START_IDX(
                x_segment_index, pic_width_in_sb, pcs_ptr->me_segments_column_count);
            x_sb_end_index = SEGMENT_END_IDX(
                x_segment_index, pic_width_in_sb, pcs_ptr->me_segments_column_count);
            y_sb_start_index = SEGMENT_START_IDX(
                y_segment_index, picture_height_in_sb, pcs_ptr->me_segments_row_count);
            y_sb_end_index = SEGMENT_END_IDX(
                y_segment_index, picture_height_in_sb, pcs_ptr->me_segments_row_count);
            // *** MOTION ESTIMATION CODE ***
            if (pcs_ptr->slice_type != I_SLICE) {
                // SB Loop
                for (y_sb_index = y_sb_start_index; y_sb_index < y_sb_end_index; ++y_sb_index) {
                    for (x_sb_index = x_sb_start_index; x_sb_index < x_sb_end_index; ++x_sb_index) {
                        sb_index    = (uint16_t)(x_sb_index + y_sb_index * pic_width_in_sb);
                        sb_origin_x = x_sb_index * scs_ptr->sb_sz;
                        sb_origin_y = y_sb_index * scs_ptr->sb_sz;

                        sb_width =
                            (pcs_ptr->aligned_width - sb_origin_x) < BLOCK_SIZE_64
                                ? pcs_ptr->aligned_width - sb_origin_x
                                : BLOCK_SIZE_64;
                        sb_height =
                            (pcs_ptr->aligned_height - sb_origin_y) < BLOCK_SIZE_64
                                ? pcs_ptr->aligned_height  - sb_origin_y
                                : BLOCK_SIZE_64;

                        // Load the SB from the input to the intermediate SB buffer
                        buffer_index = (input_picture_ptr->origin_y + sb_origin_y) *
                                           input_picture_ptr->stride_y +
                                       input_picture_ptr->origin_x + sb_origin_x;

                        context_ptr->me_context_ptr->hme_search_type = HME_RECTANGULAR;

                        for (sb_row = 0; sb_row < BLOCK_SIZE_64; sb_row++) {
                            eb_memcpy(
                                (&(context_ptr->me_context_ptr->sb_buffer[sb_row * BLOCK_SIZE_64])),
                                (&(input_picture_ptr
                                       ->buffer_y[buffer_index +
                                                  sb_row * input_picture_ptr->stride_y])),
                                BLOCK_SIZE_64 * sizeof(uint8_t));
                        }
#ifdef ARCH_X86
                        {
                            uint8_t *src_ptr = &input_padded_picture_ptr->buffer_y[buffer_index];

                            //_MM_HINT_T0     //_MM_HINT_T1    //_MM_HINT_T2//_MM_HINT_NTA
                            uint32_t i;
                            for (i = 0; i < sb_height; i++) {
                                char const *p =
                                    (char const *)(src_ptr +
                                                   i * input_padded_picture_ptr->stride_y);

                                _mm_prefetch(p, _MM_HINT_T2);

                            }
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

                            for (sb_row = 0; sb_row < (sb_height >> 1); sb_row++) {
                                eb_memcpy(
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

                            {
                                uint8_t *frame_ptr = &sixteenth_picture_ptr->buffer_y[buffer_index];
                                uint8_t *local_ptr =
                                    context_ptr->me_context_ptr->sixteenth_sb_buffer;
                                if (context_ptr->me_context_ptr->hme_search_method ==
                                    FULL_SAD_SEARCH) {
                                    for (sb_row = 0; sb_row < (sb_height >> 2); sb_row += 1) {
                                        eb_memcpy(local_ptr,
                                                  frame_ptr,
                                                  (sb_width >> 2) * sizeof(uint8_t));
                                        local_ptr += 16;
                                        frame_ptr += sixteenth_picture_ptr->stride_y;
                                    }
                                } else {
                                    for (sb_row = 0; sb_row < (sb_height >> 2); sb_row += 2) {
                                        eb_memcpy(local_ptr,
                                                  frame_ptr,
                                                  (sb_width >> 2) * sizeof(uint8_t));
                                        local_ptr += 16;
                                        frame_ptr += sixteenth_picture_ptr->stride_y << 1;
                                    }
                                }
                            }
                        }
                        context_ptr->me_context_ptr->me_alt_ref = EB_FALSE;

                        motion_estimate_sb(pcs_ptr,
                                           sb_index,
                                           sb_origin_x,
                                           sb_origin_y,
                                           context_ptr->me_context_ptr,
                                           input_picture_ptr);
                    }
                }
            }
            if (pcs_ptr->intra_pred_mode > 4)
            // *** OPEN LOOP INTRA CANDIDATE SEARCH CODE ***
            {
                // SB Loop
                for (y_sb_index = y_sb_start_index; y_sb_index < y_sb_end_index; ++y_sb_index) {
                    for (x_sb_index = x_sb_start_index; x_sb_index < x_sb_end_index; ++x_sb_index) {
                        sb_origin_x = x_sb_index * scs_ptr->sb_sz;
                        sb_origin_y = y_sb_index * scs_ptr->sb_sz;

                        sb_index = (uint16_t)(x_sb_index + y_sb_index * pic_width_in_sb);

                        open_loop_intra_search_sb(
                            pcs_ptr, sb_index, context_ptr, input_picture_ptr);
                    }
                }
            }

            // ZZ SADs Computation
            // 1 lookahead frame is needed to get valid (0,0) SAD
            if (scs_ptr->static_config.look_ahead_distance != 0) {
                // when DG is ON, the ZZ SADs are computed @ the PD process
                {
                    // ZZ SADs Computation using decimated picture
                    if (pcs_ptr->picture_number > 0) {
                        compute_decimated_zz_sad(
                            context_ptr,
                            pcs_ptr,
                            (EbPictureBufferDesc *)pa_ref_obj_
                                ->sixteenth_decimated_picture_ptr, // Hsan: always use decimated for ZZ SAD derivation until studying the trade offs and regenerating the activity threshold
                            x_sb_start_index,
                            x_sb_end_index,
                            y_sb_start_index,
                            y_sb_end_index);
                    }
                }
            }

            // Calculate the ME Distortion and OIS Historgrams

            eb_block_on_mutex(pcs_ptr->rc_distortion_histogram_mutex);

            if (scs_ptr->static_config.rate_control_mode) {
                if (pcs_ptr->slice_type != I_SLICE) {
                    uint16_t sad_interval_index;
                    for (y_sb_index = y_sb_start_index; y_sb_index < y_sb_end_index; ++y_sb_index) {
                        for (x_sb_index = x_sb_start_index; x_sb_index < x_sb_end_index;
                             ++x_sb_index) {
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

                            sb_index = (uint16_t)(x_sb_index + y_sb_index * pic_width_in_sb);
                            pcs_ptr->inter_sad_interval_index[sb_index] = 0;
                            pcs_ptr->intra_sad_interval_index[sb_index] = 0;

                            if (sb_width == BLOCK_SIZE_64 && sb_height == BLOCK_SIZE_64) {
                                sad_interval_index = (uint16_t)(
                                    pcs_ptr->rc_me_distortion[sb_index] >>
                                    (12 - SAD_PRECISION_INTERVAL)); //change 12 to 2*log2(64)

                                // SVT_LOG("%d\n", sad_interval_index);

                                sad_interval_index = (uint16_t)(sad_interval_index >> 2);
                                if (sad_interval_index > (NUMBER_OF_SAD_INTERVALS >> 1) - 1) {
                                    uint16_t sad_interval_index_temp =
                                        sad_interval_index - ((NUMBER_OF_SAD_INTERVALS >> 1) - 1);

                                    sad_interval_index = ((NUMBER_OF_SAD_INTERVALS >> 1) - 1) +
                                                         (sad_interval_index_temp >> 3);
                                }
                                if (sad_interval_index >= NUMBER_OF_SAD_INTERVALS - 1)
                                    sad_interval_index = NUMBER_OF_SAD_INTERVALS - 1;

                                pcs_ptr->inter_sad_interval_index[sb_index] = sad_interval_index;

                                pcs_ptr->me_distortion_histogram[sad_interval_index]++;

                                intra_sad_interval_index =
                                    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_64x64] >> 4;
                                intra_sad_interval_index =
                                    (uint16_t)(intra_sad_interval_index >> 2);
                                if (intra_sad_interval_index > (NUMBER_OF_SAD_INTERVALS >> 1) - 1) {
                                    uint32_t sad_interval_index_temp =
                                        intra_sad_interval_index -
                                        ((NUMBER_OF_SAD_INTERVALS >> 1) - 1);

                                    intra_sad_interval_index =
                                        ((NUMBER_OF_SAD_INTERVALS >> 1) - 1) +
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
                } else {
                    for (y_sb_index = y_sb_start_index; y_sb_index < y_sb_end_index; ++y_sb_index) {
                        for (x_sb_index = x_sb_start_index; x_sb_index < x_sb_end_index;
                             ++x_sb_index) {
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

                            sb_index = (uint16_t)(x_sb_index + y_sb_index * pic_width_in_sb);

                            pcs_ptr->inter_sad_interval_index[sb_index] = 0;
                            pcs_ptr->intra_sad_interval_index[sb_index] = 0;

                            if (sb_width == BLOCK_SIZE_64 && sb_height == BLOCK_SIZE_64) {
                                intra_sad_interval_index =
                                    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_64x64] >> 4;
                                intra_sad_interval_index =
                                    (uint16_t)(intra_sad_interval_index >> 2);
                                if (intra_sad_interval_index > (NUMBER_OF_SAD_INTERVALS >> 1) - 1) {
                                    uint32_t sad_interval_index_temp =
                                        intra_sad_interval_index -
                                        ((NUMBER_OF_SAD_INTERVALS >> 1) - 1);

                                    intra_sad_interval_index =
                                        ((NUMBER_OF_SAD_INTERVALS >> 1) - 1) +
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
                }
            }

            eb_release_mutex(pcs_ptr->rc_distortion_histogram_mutex);

            // Get Empty Results Object
            eb_get_empty_object(context_ptr->motion_estimation_results_output_fifo_ptr,
                                &out_results_wrapper_ptr);

            out_results_ptr = (MotionEstimationResults *)out_results_wrapper_ptr->object_ptr;
            out_results_ptr->pcs_wrapper_ptr = in_results_ptr->pcs_wrapper_ptr;
            out_results_ptr->segment_index   = segment_index;

            // Release the Input Results
            eb_release_object(in_results_wrapper_ptr);

            // Post the Full Results Object
            eb_post_full_object(out_results_wrapper_ptr);

        } else {
            // ME Kernel Signal(s) derivation
            tf_signal_derivation_me_kernel_oq(scs_ptr, pcs_ptr, context_ptr);

            // temporal filtering start
            context_ptr->me_context_ptr->me_alt_ref = EB_TRUE;
            svt_av1_init_temporal_filtering(
                pcs_ptr->temp_filt_pcs_list, pcs_ptr, context_ptr, in_results_ptr->segment_index);

            // Release the Input Results
            eb_release_object(in_results_wrapper_ptr);
        }
    }

    return NULL;
}
