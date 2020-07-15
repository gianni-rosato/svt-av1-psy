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
#include <emmintrin.h>
#endif
#include "EbTemporalFiltering.h"
#include "EbGlobalMotionEstimation.h"

#include "EbResize.h"

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
#if !REFACTOR_ME_HME
    uint8_t hme_me_level =
        scs_ptr->use_output_stat_file ? pcs_ptr->snd_pass_enc_mode : pcs_ptr->enc_mode;
#endif
    // HME/ME default settings
    me_context_ptr->number_hme_search_region_in_width  = 2;
    me_context_ptr->number_hme_search_region_in_height = 2;

#if !UNIFY_SC_NSC
    uint8_t sc_content_detected = pcs_ptr->sc_content_detected;
#endif

#if REFACTOR_ME_HME
    // Set the minimum ME search area
#if UNIFY_SC_NSC
    if (pcs_ptr->sc_content_detected)
#else
    if (sc_content_detected)
#endif
#if !APR02_ADOPTIONS
#if MAR30_ADOPTIONS
        if (pcs_ptr->enc_mode <= ENC_M1)
            me_context_ptr->search_area_width = me_context_ptr->search_area_height = 500;
        else
#endif
#endif
#if NEW_HME_ME_SIZES
#if PRESETS_SHIFT
#if JUNE17_ADOPTIONS
            if (pcs_ptr->enc_mode <= ENC_M3) {
#else
#if JUNE11_ADOPTIONS
            if (pcs_ptr->enc_mode <= ENC_M1) {
#else
#if MAY12_ADOPTIONS
            if (pcs_ptr->enc_mode <= ENC_M0) {
#else
#if SHIFT_M5_SC_TO_M3
            if (pcs_ptr->enc_mode <= ENC_M2) {
#else
#if APR23_ADOPTIONS_2
            if (pcs_ptr->enc_mode <= ENC_M4) {
#else
            if (pcs_ptr->enc_mode <= ENC_M2) {
#endif
#endif
#endif
#endif
#endif
#else
            if (pcs_ptr->enc_mode <= ENC_M3) {
#endif
#if NEW_M0_M1_ME_NICS
                me_context_ptr->search_area_width = me_context_ptr->search_area_height = 175;
                me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 750;
#else
                me_context_ptr->search_area_width = me_context_ptr->search_area_height = 200;
                me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 800;
#endif
            }
#if !JUNE17_ADOPTIONS
#if MAY12_ADOPTIONS
#if JUNE8_ADOPTIONS
            else if (pcs_ptr->enc_mode <= ENC_M2) {
#else
            else if (pcs_ptr->enc_mode <= ENC_M1) {
#endif
                me_context_ptr->search_area_width = me_context_ptr->search_area_height = 175;
                me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 700;
            }
#endif
#endif
#if M8_HME_ME
#if !JUNE8_ADOPTIONS
#if SHIFT_M6_SC_TO_M5
#if PRESET_SHIFITNG
            else if (pcs_ptr->enc_mode <= ENC_M2) {
#else
            else if (pcs_ptr->enc_mode <= ENC_M4) {
#endif
#else
#if APR25_3AM_ADOPTIONS
            else if (pcs_ptr->enc_mode <= ENC_M5) {
#else
#if APR24_ADOPTIONS_M6_M7
            else if (pcs_ptr->enc_mode <= ENC_M6) {
#else
            else if (pcs_ptr->enc_mode <= ENC_M5) {
#endif
#endif
#endif
#if APR23_ADOPTIONS_2
                me_context_ptr->search_area_width = me_context_ptr->search_area_height = 150;
                me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 600;
#else
                me_context_ptr->search_area_width = me_context_ptr->search_area_height = 225;
                me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 450;
#endif
            }
#endif
#if APR25_3AM_ADOPTIONS
#if APR23_4AM_M6_ADOPTIONS && !APR25_10AM_ADOPTIONS
            else if (pcs_ptr->enc_mode <= ENC_M5) {
#else
#if UNIFY_SC_NSC
            else if (pcs_ptr->enc_mode <= ENC_M6) {
#else
#if JUNE17_ADOPTIONS
            else if (pcs_ptr->enc_mode <= ENC_M5) {
#else
#if PRESET_SHIFITNG
            else if (pcs_ptr->enc_mode <= ENC_M4) {
#else
            else if (pcs_ptr->enc_mode <= ENC_M6) {
#endif
#endif
#endif
#endif
#if MAY19_ADOPTIONS
                me_context_ptr->search_area_width = me_context_ptr->search_area_height = 125;
                me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 500;
#else
                me_context_ptr->search_area_width = me_context_ptr->search_area_height = 150;
                me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 450;
#endif
            }
#endif
#if UNIFY_SC_NSC
#if NEW_M8
            else {
#else
            else if (pcs_ptr->enc_mode <= ENC_M7) {
#endif
                me_context_ptr->search_area_width = me_context_ptr->search_area_height = 75;
                me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 350;
            }
#else
#if JUNE17_ADOPTIONS
            else if (pcs_ptr->enc_mode <= ENC_M6) {
                me_context_ptr->search_area_width = me_context_ptr->search_area_height = 100;
                me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 400;
            }
#endif
#endif
#if APR23_4AM_M6_ADOPTIONS && !APR25_10AM_ADOPTIONS
            else if (pcs_ptr->enc_mode <= ENC_M6) {
                me_context_ptr->search_area_width = me_context_ptr->search_area_height = 100;
                me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 300;
            }
#endif
#if !NEW_M8
#if UPGRADE_M6_M7_M8
            else {
                me_context_ptr->search_area_width = me_context_ptr->search_area_height = 75;
                me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 300;
            }
#else
            else {
                me_context_ptr->search_area_width = me_context_ptr->search_area_height = 75;
                me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 150;
            }
#endif
#endif
#else
            else {
                me_context_ptr->search_area_width = me_context_ptr->search_area_height = 225;
                me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 450;
            }
#endif
#if NEW_M0_M1_ME_NICS
#if PRESET_SHIFITNG
    else if (pcs_ptr->enc_mode <= ENC_M1) {
#else
    else if (pcs_ptr->enc_mode <= ENC_M2) {
#endif
    me_context_ptr->search_area_width = me_context_ptr->search_area_height = 64;
#if JUNE23_ADOPTIONS
    me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 256;
#else
#if MAY17_ADOPTIONS
    me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = (input_resolution <= INPUT_SIZE_1080p_RANGE) ? 192 : 256;
#else
    me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 256;
#endif
#endif
    }
#else
#if APR22_ADOPTIONS
#if APR23_ADOPTIONS
#if M2_COMBO_1 || M1_COMBO_3 || NEW_M1_CAND
    else if (pcs_ptr->enc_mode <= ENC_M0) {
#else
    else if (pcs_ptr->enc_mode <= ENC_M2) {
#endif
#else
    else if (pcs_ptr->enc_mode <= ENC_M0) {
#endif
        me_context_ptr->search_area_width = me_context_ptr->search_area_height = 120;
        me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 480;
    }
    else if (pcs_ptr->enc_mode <= ENC_M2) {
#else
#if PRESETS_SHIFT
    else if (pcs_ptr->enc_mode <= ENC_M2) {
#else
    else if (pcs_ptr->enc_mode <= ENC_M3) {
#endif
#endif
        me_context_ptr->search_area_width = me_context_ptr->search_area_height = 80;
        me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 320;
    }
#endif
#if M8_HME_ME
#if JUNE23_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M2) {
        me_context_ptr->search_area_width = me_context_ptr->search_area_height = 64;
        me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 192;
    }
#else
#if JUNE8_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M1) {
#else
#if NEW_M5_HME_ME
#if PRESET_SHIFITNG
    else if (pcs_ptr->enc_mode <= ENC_M2) {
#else
    else if (pcs_ptr->enc_mode <= ENC_M4) {
#endif
#else
    else if (pcs_ptr->enc_mode <= ENC_M5) {
#endif
#endif
#if MAY12_ADOPTIONS
    me_context_ptr->search_area_width = me_context_ptr->search_area_height = 64;
#if MAY17_ADOPTIONS
    me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = (input_resolution <= INPUT_SIZE_1080p_RANGE) ? 192 : 256;
#else
    me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 256;
#endif
#else
        me_context_ptr->search_area_width = me_context_ptr->search_area_height = 75;
        me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 150;
#endif
    }
#endif
#if NEW_M5_HME_ME
#if PRESET_SHIFITNG
    else if (pcs_ptr->enc_mode <= ENC_M3) {
#else
    else if (pcs_ptr->enc_mode <= ENC_M5) {
#endif
        me_context_ptr->search_area_width = me_context_ptr->search_area_height = 64;
        me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 128;
    }
#endif
#if APR25_1PM_ADOPTIONS
#if JUNE25_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M4) {
#else
#if JUNE17_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M5) {
#else
#if PRESET_SHIFITNG
    else if (pcs_ptr->enc_mode <= ENC_M4) {
#else
    else if (pcs_ptr->enc_mode <= ENC_M6) {
#endif
#endif
#endif
    me_context_ptr->search_area_width = me_context_ptr->search_area_height = 32;
    me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 128;
    }
#endif
#if UPGRADE_M6_M7_M8
#if !APR25_11AM_ADOPTIONS
#if APR25_10AM_ADOPTIONS
    else if (pcs_ptr->enc_mode <= ENC_M6) {
#else
    else if (pcs_ptr->enc_mode <= ENC_M7) {
#endif
        me_context_ptr->search_area_width = me_context_ptr->search_area_height = 25;
        me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 100;
    }
#endif
    else {
        me_context_ptr->search_area_width = me_context_ptr->search_area_height = 16;
        me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 64;
    }
#else
    else {
        me_context_ptr->search_area_width = me_context_ptr->search_area_height = 25;
        me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 50;
    }
#endif
#else
    else {
        me_context_ptr->search_area_width = me_context_ptr->search_area_height = 75;
        me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = 150;
    }
#endif
#else
        if (pcs_ptr->enc_mode <= ENC_M3)
            me_context_ptr->search_area_width = me_context_ptr->search_area_height = 390;
        else
            me_context_ptr->search_area_width = me_context_ptr->search_area_height = 225;
    else if (pcs_ptr->enc_mode <= ENC_M3)
#if MAR30_ADOPTIONS
        me_context_ptr->search_area_width = me_context_ptr->search_area_height = 120;
#else
        me_context_ptr->search_area_width = me_context_ptr->search_area_height = 150;
#endif
    else
        me_context_ptr->search_area_width = me_context_ptr->search_area_height = 75;

    me_context_ptr->max_me_search_width = me_context_ptr->search_area_width * 2;
    me_context_ptr->max_me_search_height = me_context_ptr->search_area_height * 2;

    me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = me_context_ptr->max_me_search_height;
#endif

#if ADD_MAX_HME_SIGNAL
#if NEW_HME_ME_SIZES
#if !UNIFY_SC_NSC
    if (sc_content_detected) {
#if UPGRADE_M6_M7_M8
#if MAY19_ADOPTIONS
#if PRESET_SHIFITNG
        if (pcs_ptr->enc_mode <= ENC_M3) {
#else
        if (pcs_ptr->enc_mode <= ENC_M5) {
#endif
#else
#if APR25_10AM_ADOPTIONS
        if (pcs_ptr->enc_mode <= ENC_M6) {
#else
        if (pcs_ptr->enc_mode <= ENC_M7) {
#endif
#endif
#if MAY16_M0_ADOPTIONS
            me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = 240;
            me_context_ptr->hme_level0_max_total_search_area_width = me_context_ptr->hme_level0_max_total_search_area_height = 480;
#else
#if HME_4K_ADOPTIONS
            me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = input_resolution <= INPUT_SIZE_1080p_RANGE ? 120 : 240;
            me_context_ptr->hme_level0_max_total_search_area_width = me_context_ptr->hme_level0_max_total_search_area_height = 480;
#else
            me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = 150;
            me_context_ptr->hme_level0_max_total_search_area_width = me_context_ptr->hme_level0_max_total_search_area_height = 600;
#endif
#endif
        }
        else {
            me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = 75;
            me_context_ptr->hme_level0_max_total_search_area_width = me_context_ptr->hme_level0_max_total_search_area_height = 300;
        }
#else
        me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = 150;
        me_context_ptr->hme_level0_max_total_search_area_width = me_context_ptr->hme_level0_max_total_search_area_height = 600;
#endif
    }
    else {
#endif
#if APR22_ADOPTIONS
#if MAY19_ADOPTIONS
#if JUNE23_ADOPTIONS
        if (pcs_ptr->enc_mode <= ENC_M4) {
#else
#if PRESET_SHIFITNG
        if (pcs_ptr->enc_mode <= ENC_M3) {
#else
        if (pcs_ptr->enc_mode <= ENC_M5) {
#endif
#endif
#else
#if MAY12_ADOPTIONS
        if (pcs_ptr->enc_mode <= ENC_M4) {
#else
#if APR23_ADOPTIONS
#if M1_COMBO_3 || M2_COMBO_2 || M2_COMBO_3 || NEW_M1_CAND
        if (pcs_ptr->enc_mode <= ENC_M0) {
#else
        if (pcs_ptr->enc_mode <= ENC_M2) {
#endif
#else
        if (pcs_ptr->enc_mode <= ENC_M0) {
#endif
#endif
#endif
#if HME_4K_ADOPTIONS
            me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = input_resolution <= INPUT_SIZE_1080p_RANGE ? 120 : 240;
            me_context_ptr->hme_level0_max_total_search_area_width = me_context_ptr->hme_level0_max_total_search_area_height = 480;
#else
            me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = 120;
            me_context_ptr->hme_level0_max_total_search_area_width = me_context_ptr->hme_level0_max_total_search_area_height = 480;
#endif
        }
#if UNIFY_SC_NSC
#if !JUNE23_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M4) {
            me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = 75;
            me_context_ptr->hme_level0_max_total_search_area_width = me_context_ptr->hme_level0_max_total_search_area_height = 300;
        }
#endif
#if !JUNE25_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M5) {
            me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = 48;
            me_context_ptr->hme_level0_max_total_search_area_width = me_context_ptr->hme_level0_max_total_search_area_height = 192;
        }
#endif
        else if (pcs_ptr->enc_mode <= ENC_M6) {
            me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = 32;
            me_context_ptr->hme_level0_max_total_search_area_width = me_context_ptr->hme_level0_max_total_search_area_height = 164;
        }
#endif
#if UPGRADE_M6_M7_M8
#if !MAY12_ADOPTIONS
#if NEW_M5_HME_ME
        else if (pcs_ptr->enc_mode <= ENC_M4) {
#else
#if APR25_1PM_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M5) {
#else
#if APR25_12AM_ADOPTIONS
        else if (pcs_ptr->enc_mode <= ENC_M6) {
#else
        else if (pcs_ptr->enc_mode <= ENC_M7) {
#endif
#endif
#endif
            me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = 100;
            me_context_ptr->hme_level0_max_total_search_area_width = me_context_ptr->hme_level0_max_total_search_area_height = 400;
        }
#endif
#if !MAY19_ADOPTIONS
#if NEW_M5_HME_ME
        else if (pcs_ptr->enc_mode <= ENC_M5) {
            me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = 64;
            me_context_ptr->hme_level0_max_total_search_area_width = me_context_ptr->hme_level0_max_total_search_area_height = 256;
        }
#endif
#endif
        else {
            me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = 32;
            me_context_ptr->hme_level0_max_total_search_area_width = me_context_ptr->hme_level0_max_total_search_area_height = 128;
        }
#else
        else {
            me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = 100;
            me_context_ptr->hme_level0_max_total_search_area_width = me_context_ptr->hme_level0_max_total_search_area_height = 400;
        }
#endif
#else
        me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = 100;
        me_context_ptr->hme_level0_max_total_search_area_width = me_context_ptr->hme_level0_max_total_search_area_height = 400;
#endif
#if !UNIFY_SC_NSC
    }
#endif
#else
    me_context_ptr->hme_level0_max_total_search_area_width =
        me_context_ptr->hme_level0_max_total_search_area_height =
        me_context_ptr->hme_level0_total_search_area_width * 3;
#endif

    me_context_ptr->hme_level0_max_search_area_in_width_array[0] =
        me_context_ptr->hme_level0_max_search_area_in_width_array[1] =
        me_context_ptr->hme_level0_max_total_search_area_width / me_context_ptr->number_hme_search_region_in_width;
    me_context_ptr->hme_level0_max_search_area_in_height_array[0] =
        me_context_ptr->hme_level0_max_search_area_in_height_array[1] =
        me_context_ptr->hme_level0_max_total_search_area_height / me_context_ptr->number_hme_search_region_in_height;
#endif
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

#if ADD_HME_DECIMATION_SIGNAL
#if MAR30_ADOPTIONS
#if !UNIFY_SC_NSC
    if (pcs_ptr->sc_content_detected)
#if APR02_ADOPTIONS
        me_context_ptr->hme_decimation = TWO_DECIMATION_HME;
#else
        me_context_ptr->hme_decimation = pcs_ptr->enc_mode <= ENC_M1 ? ONE_DECIMATION_HME : TWO_DECIMATION_HME;
#endif
    else if (input_resolution <= INPUT_SIZE_720p_RANGE)
#else
    if (input_resolution <= INPUT_SIZE_720p_RANGE)
#endif
#if PRESETS_SHIFT
#if PRESET_SHIFITNG
        me_context_ptr->hme_decimation = pcs_ptr->enc_mode <= ENC_M1 ? ONE_DECIMATION_HME : TWO_DECIMATION_HME;
#else
        me_context_ptr->hme_decimation = pcs_ptr->enc_mode <= ENC_M2 ? ONE_DECIMATION_HME : TWO_DECIMATION_HME;
#endif
#else
        me_context_ptr->hme_decimation = pcs_ptr->enc_mode <= ENC_M3 ? ONE_DECIMATION_HME : TWO_DECIMATION_HME;
#endif
    else
        me_context_ptr->hme_decimation = TWO_DECIMATION_HME;
#else
    me_context_ptr->hme_decimation = TWO_DECIMATION_HME;
#endif
#endif

    // Scale up the MIN ME area if low frame rate
#if UNIFY_SC_NSC
    uint8_t  low_frame_rate_flag = (scs_ptr->static_config.frame_rate >> 16) < 50 ? 1 : 0;
#else
    uint8_t  low_frame_rate_flag = sc_content_detected ? 0 :
        (scs_ptr->static_config.frame_rate >> 16) < 50 ? 1 : 0;
#endif
    if (low_frame_rate_flag) {
        me_context_ptr->search_area_width = (me_context_ptr->search_area_width * 3) / 2;
        me_context_ptr->search_area_height = (me_context_ptr->search_area_height * 3) / 2;
    }
#else
    // sequence frame rate based mutiplier
    uint8_t  fr_rate_mult_num   = sc_content_detected ? 1 :
        (scs_ptr->static_config.frame_rate >> 16) < 50 ? 3 : 1;
    uint8_t  fr_rate_mult_denum = sc_content_detected ? 1 :
        (scs_ptr->static_config.frame_rate >> 16) < 50 ? 2 : 1;

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

    me_context_ptr->search_area_width  =
        (min_me_search_width[sc_content_detected][input_resolution][hme_me_level]
            *fr_rate_mult_num)/fr_rate_mult_denum ;
    me_context_ptr->search_area_height =
        (min_me_search_height[sc_content_detected][input_resolution][hme_me_level]
            *fr_rate_mult_num)/fr_rate_mult_denum ;

#endif

#if !REMOVE_ME_SUBPEL_CODE
    assert(me_context_ptr->search_area_width <= MAX_SEARCH_AREA_WIDTH &&
           "increase MAX_SEARCH_AREA_WIDTH");
    assert(me_context_ptr->search_area_height <= MAX_SEARCH_AREA_HEIGHT &&
           "increase MAX_SEARCH_AREA_HEIGHT");
#if ME_MEM_OPT2
    assert(me_context_ptr->max_me_search_width <= MAX_SEARCH_AREA_WIDTH && "increase MAX_SEARCH_AREA_WIDTH");
    assert(me_context_ptr->max_me_search_height <= MAX_SEARCH_AREA_HEIGHT && "increase MAX_SEARCH_AREA_HEIGHT");
#endif
#endif

    me_context_ptr->update_hme_search_center_flag = 1;

#if NEW_RESOLUTION_RANGES
    if (input_resolution <= INPUT_SIZE_480p_RANGE)
#else
    if (input_resolution <= INPUT_SIZE_576p_RANGE_OR_LOWER)
#endif
        me_context_ptr->update_hme_search_center_flag = 0;

    return NULL;
};
#if ME_HME_PRUNING_CLEANUP
void set_me_hme_ref_prune_ctrls(MeContext* context_ptr, uint8_t prune_level) {
    MeHmeRefPruneCtrls* me_hme_prune_ctrls = &context_ptr->me_hme_prune_ctrls;

    switch (prune_level)
    {
    case 0:
#if ON_OFF_FEATURE_MRP
        me_hme_prune_ctrls->enable_me_hme_ref_pruning = override_feature_level (context_ptr->mrp_level, 0, 0, 1);
#else
        me_hme_prune_ctrls->enable_me_hme_ref_pruning = 0;
#endif
        me_hme_prune_ctrls->prune_ref_if_hme_sad_dev_bigger_than_th = (uint16_t)~0;
        me_hme_prune_ctrls->prune_ref_if_me_sad_dev_bigger_than_th = (uint16_t)~0;
        break;
    case 1:
#if ON_OFF_FEATURE_MRP
        me_hme_prune_ctrls->enable_me_hme_ref_pruning = override_feature_level(context_ptr->mrp_level, 1, 0, 1);
#else
        me_hme_prune_ctrls->enable_me_hme_ref_pruning = 1;
#endif
        me_hme_prune_ctrls->prune_ref_if_hme_sad_dev_bigger_than_th = 160;
        me_hme_prune_ctrls->prune_ref_if_me_sad_dev_bigger_than_th = (uint16_t)~0;
        break;
    case 2:
#if ON_OFF_FEATURE_MRP
        me_hme_prune_ctrls->enable_me_hme_ref_pruning = override_feature_level(context_ptr->mrp_level, 1, 0, 1);
#else
        me_hme_prune_ctrls->enable_me_hme_ref_pruning = 1;
#endif
        me_hme_prune_ctrls->prune_ref_if_hme_sad_dev_bigger_than_th = 80;
        me_hme_prune_ctrls->prune_ref_if_me_sad_dev_bigger_than_th = 60;
        break;
    case 3:
#if ON_OFF_FEATURE_MRP
        me_hme_prune_ctrls->enable_me_hme_ref_pruning = override_feature_level(context_ptr->mrp_level, 1, 0, 1);
#else
        me_hme_prune_ctrls->enable_me_hme_ref_pruning = 1;
#endif
        me_hme_prune_ctrls->prune_ref_if_hme_sad_dev_bigger_than_th = 50;
        me_hme_prune_ctrls->prune_ref_if_me_sad_dev_bigger_than_th = 60;
        break;
    case 4:
#if ON_OFF_FEATURE_MRP
        me_hme_prune_ctrls->enable_me_hme_ref_pruning = override_feature_level(context_ptr->mrp_level, 1, 0, 1);
#else
        me_hme_prune_ctrls->enable_me_hme_ref_pruning = 1;
#endif
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
#endif
/******************************************************
* Derive ME Settings for OQ
  Input   : encoder mode and tune
  Output  : ME Kernel signal(s)
******************************************************/
EbErrorType signal_derivation_me_kernel_oq(SequenceControlSet *       scs_ptr,
                                           PictureParentControlSet *  pcs_ptr,
                                           MotionEstimationContext_t *context_ptr) {
    EbErrorType return_error = EB_ErrorNone;

#if REMOVE_MR_MACRO
    EbEncMode enc_mode = scs_ptr->use_output_stat_file ?
        pcs_ptr->snd_pass_enc_mode : pcs_ptr->enc_mode;
#else
    uint8_t enc_mode =
        scs_ptr->use_output_stat_file ? pcs_ptr->snd_pass_enc_mode : pcs_ptr->enc_mode;
#endif
#if ON_OFF_FEATURE_MRP
    context_ptr->me_context_ptr->mrp_level = pcs_ptr->mrp_level;
#endif
    EbInputResolution input_resolution = scs_ptr->input_resolution;
    uint8_t sc_content_detected = pcs_ptr->sc_content_detected;
#if !REFACTOR_ME_HME
    uint8_t  hme_me_level = scs_ptr->use_output_stat_file ?
        pcs_ptr->snd_pass_enc_mode : pcs_ptr->enc_mode;
#endif
    // Set ME/HME search regions
    if (scs_ptr->static_config.use_default_me_hme)
        set_me_hme_params_oq(
            context_ptr->me_context_ptr, pcs_ptr, scs_ptr, input_resolution);

    else
        set_me_hme_params_from_config(scs_ptr, context_ptr->me_context_ptr);

#if !REFACTOR_ME_HME
    // ME
    context_ptr->me_context_ptr->max_me_search_width =
        max_me_search_width[sc_content_detected][input_resolution][hme_me_level];
    context_ptr->me_context_ptr->max_me_search_height =
        max_me_search_height[sc_content_detected][input_resolution][hme_me_level];
#endif
#if !REMOVE_ME_SUBPEL_CODE
    if (sc_content_detected)
#if MAR11_ADOPTIONS
        // fractional_search_method is not used if subpel is OFF
        context_ptr->me_context_ptr->fractional_search_method = FULL_SAD_SEARCH;
#else
        context_ptr->me_context_ptr->fractional_search_method =
            (enc_mode <= ENC_M1) ? FULL_SAD_SEARCH : SUB_SAD_SEARCH;
#endif
    else
#if MAR2_M8_ADOPTIONS
        context_ptr->me_context_ptr->fractional_search_method = SSD_SEARCH;
#else
        if (enc_mode <= ENC_M6)
        context_ptr->me_context_ptr->fractional_search_method = SSD_SEARCH;
    else
        context_ptr->me_context_ptr->fractional_search_method = FULL_SAD_SEARCH;
#endif
#endif
    // Set HME flags
    context_ptr->me_context_ptr->enable_hme_flag        = pcs_ptr->enable_hme_flag;
    context_ptr->me_context_ptr->enable_hme_level0_flag = pcs_ptr->enable_hme_level0_flag;
    context_ptr->me_context_ptr->enable_hme_level1_flag = pcs_ptr->enable_hme_level1_flag;
    context_ptr->me_context_ptr->enable_hme_level2_flag = pcs_ptr->enable_hme_level2_flag;
#if !REMOVE_ME_SUBPEL_CODE
    if (scs_ptr->static_config.enable_subpel == DEFAULT)
        // Set the default settings of subpel
        if (sc_content_detected)
                context_ptr->me_context_ptr->use_subpel_flag = 0;
        else
            context_ptr->me_context_ptr->use_subpel_flag = 1;
    else
        context_ptr->me_context_ptr->use_subpel_flag = scs_ptr->static_config.enable_subpel;
#if M1_C2_ADOPTIONS
    if (enc_mode <= ENC_M0) {
#else
#if APR23_ADOPTIONS
    if (enc_mode <= ENC_M2) {
#else
#if MAR10_ADOPTIONS
    if (enc_mode <= ENC_M1) {
#else
    if (enc_mode <= ENC_M0) {
#endif
#endif
#endif
        context_ptr->me_context_ptr->half_pel_mode =
            (sc_content_detected) ? REFINEMENT_HP_MODE : EX_HP_MODE;
    }

#if !MAR11_ADOPTIONS
#if MAR3_M2_ADOPTIONS
    else if (enc_mode <= ENC_M2) {
#else
    else if (enc_mode <= ENC_M2) {
#endif
        context_ptr->me_context_ptr->half_pel_mode =
            (sc_content_detected) ? REFINEMENT_HP_MODE : SWITCHABLE_HP_MODE;
    }
#endif
    else {
        context_ptr->me_context_ptr->half_pel_mode    = REFINEMENT_HP_MODE;
    }
    context_ptr->me_context_ptr->h_pel_search_wind = H_PEL_SEARCH_WIND_2;
    // Set fractional search model
    // 0: search all blocks
    // 1: selective based on Full-Search SAD & MV.
    // 2: off
    if (context_ptr->me_context_ptr->use_subpel_flag == 1) {
#if MAR17_ADOPTIONS
        context_ptr->me_context_ptr->fractional_search_model = 0;
#else
        if (enc_mode <= ENC_M6)
            context_ptr->me_context_ptr->fractional_search_model = 0;
        else
            context_ptr->me_context_ptr->fractional_search_model = 1;
#endif
    } else
        context_ptr->me_context_ptr->fractional_search_model = 2;
#endif
    // HME Search Method
    context_ptr->me_context_ptr->hme_search_method = SUB_SAD_SEARCH;
    // ME Search Method
    context_ptr->me_context_ptr->me_search_method = SUB_SAD_SEARCH;

    if (scs_ptr->static_config.enable_global_motion == EB_TRUE &&
        pcs_ptr->frame_superres_enabled == EB_FALSE) {
#if UNIFY_SC_NSC
#if JUNE26_ADOPTIONS
        if (enc_mode <= ENC_M6)
#else
        if (enc_mode <= ENC_M5)
#endif
            context_ptr->me_context_ptr->compute_global_motion = EB_TRUE;
        else
            context_ptr->me_context_ptr->compute_global_motion = EB_FALSE;
#else
#if MAR4_M6_ADOPTIONS
        if (pcs_ptr->sc_content_detected)
#if MAY19_ADOPTIONS
#if JUNE17_ADOPTIONS
            if (enc_mode <= ENC_M5)
#else
#if PRESET_SHIFITNG
            if (enc_mode <= ENC_M4)
#else
            if (enc_mode <= ENC_M6)
#endif
#endif
#else
#if MAY12_ADOPTIONS
            if (enc_mode <= ENC_M4)
#else
#if SHIFT_M5_SC_TO_M3
            if (enc_mode <= ENC_M2)
#else
#if PRESETS_SHIFT
            if (enc_mode <= ENC_M4)
#else
#if MAR12_ADOPTIONS
#if MAR17_ADOPTIONS
            if (enc_mode <= ENC_M7)
#else
            if (enc_mode <= ENC_M3)
#endif
#else
#if MAR10_ADOPTIONS
#if MAR11_ADOPTIONS
            if (enc_mode <= ENC_M1)
#else
            if (enc_mode <= ENC_M2)
#endif
#else
            if (enc_mode <= ENC_M3)
#endif
#endif
#endif
#endif
#endif
#endif
                context_ptr->me_context_ptr->compute_global_motion = EB_TRUE;
            else
                context_ptr->me_context_ptr->compute_global_motion = EB_FALSE;
#if MAY19_ADOPTIONS
#if JUNE17_ADOPTIONS
        else if (enc_mode <= ENC_M5)
#else
#if PRESET_SHIFITNG
        else if (enc_mode <= ENC_M4)
#else
        else if (enc_mode <= ENC_M6)
#endif
#endif
#else
#if PRESETS_SHIFT
        else if (enc_mode <= ENC_M4)
#else
#if MAR17_ADOPTIONS
        else if (enc_mode <= ENC_M7)
#else
        else if (enc_mode <= ENC_M5)
#endif
#endif
#endif
#else
        if (enc_mode <= ENC_M1)
#endif
            context_ptr->me_context_ptr->compute_global_motion = EB_TRUE;
        else
            context_ptr->me_context_ptr->compute_global_motion = EB_FALSE;
#endif
#if GM_LIST1
        //TODO: enclose all gm signals into a control
        context_ptr->me_context_ptr->gm_identiy_exit = EB_FALSE;
#endif
    } else
        context_ptr->me_context_ptr->compute_global_motion = EB_FALSE;

#if !SHUT_ME_NSQ_SEARCH
    // Me nsq search levels.
    // 0: feature off -> perform nsq_search.
    // 1: perform me nsq_search only for the best refrenece picture.
    // 2: perform me nsq_search only for the nearest refrenece pictures.
    // 3: me nsq_search off.
#if APR23_ADOPTIONS
    if (MR_MODE)
#else
#if MAR30_ADOPTIONS
    if (enc_mode <= ENC_M0)
#else
#if MAR10_ADOPTIONS
    if (enc_mode <= ENC_M1 && pcs_ptr->sc_content_detected == 0)
#else
    if (MR_MODE && pcs_ptr->sc_content_detected == 0)
#endif
#endif
#endif
        context_ptr->me_context_ptr->inherit_rec_mv_from_sq_block = 0;
#if UPGRADE_M8
#if MAY19_ADOPTIONS
    else if (enc_mode <= ENC_M5 || (!pcs_ptr->sc_content_detected && enc_mode <= ENC_M6))
#else
#if NSQ_OFF_IN_M6_M7_ME
    else if (enc_mode <= ENC_M5)
#else
    else if (enc_mode <= ENC_M7)
#endif
#endif
        context_ptr->me_context_ptr->inherit_rec_mv_from_sq_block = 2;
    else
        context_ptr->me_context_ptr->inherit_rec_mv_from_sq_block = 3;
#else
    else
        context_ptr->me_context_ptr->inherit_rec_mv_from_sq_block = 2;
#endif
#endif
#if ME_HME_PRUNING_CLEANUP

    // Set hme/me based reference pruning level (0-4)
#if M0_HME_ME_PRUNE
#if !UNIFY_SC_NSC
    if (pcs_ptr->sc_content_detected)
        if (enc_mode <= ENC_M0)
            set_me_hme_ref_prune_ctrls(context_ptr->me_context_ptr, 0);
        else
            set_me_hme_ref_prune_ctrls(context_ptr->me_context_ptr, 1);
    else
#endif
#if JUNE23_ADOPTIONS
#if REMOVE_MR_MACRO
        if (enc_mode <= ENC_MR)
#else
        if (MR_MODE)
#endif
#else
        if (enc_mode <= ENC_M0)
#endif
            set_me_hme_ref_prune_ctrls(context_ptr->me_context_ptr, 0);
#if JUNE11_ADOPTIONS
        else if (enc_mode <= ENC_M3)
#else
        else if (enc_mode <= ENC_M2)
#endif
            set_me_hme_ref_prune_ctrls(context_ptr->me_context_ptr, 2);
        else
            set_me_hme_ref_prune_ctrls(context_ptr->me_context_ptr, 4);
#else
    if (MR_MODE)
        set_me_hme_ref_prune_ctrls(context_ptr->me_context_ptr, 0);
    else if (pcs_ptr->sc_content_detected)
#if NEW_MRP_SETTINGS
        if (enc_mode <= ENC_M0)
            set_me_hme_ref_prune_ctrls(context_ptr->me_context_ptr, 0);
        else
            set_me_hme_ref_prune_ctrls(context_ptr->me_context_ptr, 1);
#else
#if ADOPT_SC_HME_PRUNING
        set_me_hme_ref_prune_ctrls(context_ptr->me_context_ptr, 1);
#else
        set_me_hme_ref_prune_ctrls(context_ptr->me_context_ptr, 0);
#endif
#endif
#if JUNE8_ADOPTIONS
    else if (enc_mode <= ENC_M2)
#else
#if MAY12_ADOPTIONS
#if PRESET_SHIFITNG
    else if (enc_mode <= ENC_M1)
#else
    else if (enc_mode <= ENC_M2)
#endif
#else
#if M1_COMBO_1 || NEW_M1_CAND
    else if (enc_mode <= ENC_M0)
#else
    else if (enc_mode <= ENC_M1)
#endif
#endif
#endif
        set_me_hme_ref_prune_ctrls(context_ptr->me_context_ptr, 2);
#if !MAY12_ADOPTIONS
    else if (enc_mode <= ENC_M2)
        set_me_hme_ref_prune_ctrls(context_ptr->me_context_ptr, 3);
#endif
    else
        set_me_hme_ref_prune_ctrls(context_ptr->me_context_ptr, 4);
#endif
    // Set hme-based me sr adjustment level
#if !MAY19_ADOPTIONS
    if (MR_MODE)
        set_me_sr_adjustment_ctrls(context_ptr->me_context_ptr, 0);
    else if (pcs_ptr->sc_content_detected)
#else
#if JUNE15_ADOPTIONS
#if REMOVE_MR_MACRO
    if (enc_mode <= ENC_MRS)
#else
    if (MRS_MODE)
#endif
        set_me_sr_adjustment_ctrls(context_ptr->me_context_ptr, 0);
#if !UNIFY_SC_NSC
    else if (pcs_ptr->sc_content_detected)
#endif
#else
    if (pcs_ptr->sc_content_detected)
#endif
#endif
#if !UNIFY_SC_NSC
        set_me_sr_adjustment_ctrls(context_ptr->me_context_ptr, 1);
#endif
    else
        set_me_sr_adjustment_ctrls(context_ptr->me_context_ptr, 2);
#else
#if ADD_ME_SIGNAL_FOR_PRUNING_TH
#if MAR20_ADOPTIONS
    if (pcs_ptr->sc_content_detected)
#if MAR30_ADOPTIONS
        context_ptr->me_context_ptr->prune_ref_if_hme_sad_dev_bigger_than_th = 300;
#else
        context_ptr->me_context_ptr->prune_ref_if_hme_sad_dev_bigger_than_th = 200;
#endif
    else if (enc_mode <= ENC_M1)
        context_ptr->me_context_ptr->prune_ref_if_hme_sad_dev_bigger_than_th = 80;
#if PRESETS_SHIFT
    else if (enc_mode <= ENC_M2)
#else
    else if (enc_mode <= ENC_M3)
#endif
        context_ptr->me_context_ptr->prune_ref_if_hme_sad_dev_bigger_than_th = 50;
#else
#if MAR19_ADOPTIONS
    if (enc_mode <= ENC_M4)
#else
    if (enc_mode <= ENC_M7)
#endif
        context_ptr->me_context_ptr->prune_ref_if_hme_sad_dev_bigger_than_th = 80;
#endif
    else
        context_ptr->me_context_ptr->prune_ref_if_hme_sad_dev_bigger_than_th = 30;

#if APR08_ADOPTIONS
    if (MR_MODE)
        context_ptr->me_context_ptr->prune_ref_if_me_sad_dev_bigger_than_th = (uint16_t)~0;
    else
        context_ptr->me_context_ptr->prune_ref_if_me_sad_dev_bigger_than_th = 60;
#else
#if MAR19_ADOPTIONS
    // Turn this threshold off
    context_ptr->me_context_ptr->prune_ref_if_me_sad_dev_bigger_than_th = (uint16_t)~0;
#else
    if (enc_mode <= ENC_M7)
        context_ptr->me_context_ptr->prune_ref_if_me_sad_dev_bigger_than_th = 30;
    else
        context_ptr->me_context_ptr->prune_ref_if_me_sad_dev_bigger_than_th = 15;
#endif
#endif
#endif
#endif
#if !ADD_MAX_HME_SIGNAL
#if ADD_HME_MIN_MAX_MULTIPLIER_SIGNAL
    context_ptr->me_context_ptr->max_hme_sr_area_multipler = 3;
#endif
#endif
    return return_error;
};

/************************************************
 * Set ME/HME Params for Altref Temporal Filtering
 ************************************************/
void *tf_set_me_hme_params_oq(MeContext *me_context_ptr, PictureParentControlSet *pcs_ptr,
                              SequenceControlSet *scs_ptr, EbInputResolution input_resolution) {
    UNUSED(scs_ptr);
#if REFACTOR_ME_HME && MAR12_ADOPTIONS
    UNUSED(pcs_ptr);
#endif
#if !REFACTOR_ME_HME || !MAR12_ADOPTIONS
    uint8_t hme_me_level =
        scs_ptr->use_output_stat_file ? pcs_ptr->snd_pass_enc_mode : pcs_ptr->enc_mode;
#endif
    // HME/ME default settings
    me_context_ptr->number_hme_search_region_in_width  = 2;
    me_context_ptr->number_hme_search_region_in_height = 2;
#if !REFACTOR_ME_HME
    uint8_t sc_content_detected = pcs_ptr->sc_content_detected;
#endif
#if REFACTOR_ME_HME
    // Set the minimum ME search area
#if APR22_ADOPTIONS
    me_context_ptr->search_area_width = me_context_ptr->search_area_height = input_resolution <= INPUT_SIZE_480p_RANGE ? 60 : 16;
    me_context_ptr->max_me_search_width = me_context_ptr->max_me_search_height = input_resolution <= INPUT_SIZE_480p_RANGE ? 120 : 32;
#else
    me_context_ptr->search_area_width = me_context_ptr->search_area_height = 16;

    me_context_ptr->max_me_search_width = me_context_ptr->search_area_width * 2;
    me_context_ptr->max_me_search_height = me_context_ptr->search_area_height * 2;
#endif
#if APR22_ADOPTIONS
#if UPGRADE_M6_M7_M8
#if JUNE26_ADOPTIONS
    if (pcs_ptr->enc_mode <= ENC_M5)
#else
#if JUNE25_ADOPTIONS
    if (pcs_ptr->enc_mode <= ENC_M4)
#else
#if PRESET_SHIFITNG
    if (pcs_ptr->enc_mode <= ENC_M5)
#else
    if (pcs_ptr->enc_mode <= ENC_M7)
#endif
#endif
#endif
        me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = 30;
    else
        me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = 16;
#else
    me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = 30;
#endif
#else
    me_context_ptr->hme_level0_total_search_area_width = me_context_ptr->hme_level0_total_search_area_height = me_context_ptr->max_me_search_height;
#endif
#if ADD_MAX_HME_SIGNAL
#if APR22_ADOPTIONS
#if UPGRADE_M6_M7_M8
#if JUNE26_ADOPTIONS
    if (pcs_ptr->enc_mode <= ENC_M5)
#else
#if JUNE25_ADOPTIONS
    if (pcs_ptr->enc_mode <= ENC_M4)
#else
#if PRESET_SHIFITNG
    if (pcs_ptr->enc_mode <= ENC_M5)
#else
    if (pcs_ptr->enc_mode <= ENC_M7)
#endif
#endif
#endif
        me_context_ptr->hme_level0_max_total_search_area_width = me_context_ptr->hme_level0_max_total_search_area_height = 60;
    else
        me_context_ptr->hme_level0_max_total_search_area_width = me_context_ptr->hme_level0_max_total_search_area_height = 32;
#else
    me_context_ptr->hme_level0_max_total_search_area_width =
        me_context_ptr->hme_level0_max_total_search_area_height = 60;
#endif
#else
    me_context_ptr->hme_level0_max_total_search_area_width =
        me_context_ptr->hme_level0_max_total_search_area_height =
        me_context_ptr->hme_level0_total_search_area_width * 3;
#endif


    me_context_ptr->hme_level0_max_search_area_in_width_array[0] =
        me_context_ptr->hme_level0_max_search_area_in_width_array[1] =
        me_context_ptr->hme_level0_max_total_search_area_width / me_context_ptr->number_hme_search_region_in_width;
    me_context_ptr->hme_level0_max_search_area_in_height_array[0] =
        me_context_ptr->hme_level0_max_search_area_in_height_array[1] =
        me_context_ptr->hme_level0_max_total_search_area_height / me_context_ptr->number_hme_search_region_in_height;
#endif
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

#if ADD_HME_DECIMATION_SIGNAL
    me_context_ptr->hme_decimation = TWO_DECIMATION_HME;
#endif
#else
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
    me_context_ptr->search_area_width  =
        min_metf_search_width[sc_content_detected][input_resolution][hme_me_level];
    me_context_ptr->search_area_height =
        min_metf_search_height[sc_content_detected][input_resolution][hme_me_level];
#endif
#if !REMOVE_ME_SUBPEL_CODE
    assert(me_context_ptr->search_area_width <= MAX_SEARCH_AREA_WIDTH &&
           "increase MAX_SEARCH_AREA_WIDTH");
    assert(me_context_ptr->search_area_height <= MAX_SEARCH_AREA_HEIGHT &&
           "increase MAX_SEARCH_AREA_HEIGHT");
#if ME_MEM_OPT2
    assert(me_context_ptr->max_me_search_width <= MAX_SEARCH_AREA_WIDTH && "increase MAX_SEARCH_AREA_WIDTH");
    assert(me_context_ptr->max_me_search_height <= MAX_SEARCH_AREA_HEIGHT && "increase MAX_SEARCH_AREA_HEIGHT");
#endif
#endif
    me_context_ptr->update_hme_search_center_flag = 1;

#if NEW_RESOLUTION_RANGES
    if (input_resolution <= INPUT_SIZE_480p_RANGE)
#else
    if (input_resolution <= INPUT_SIZE_576p_RANGE_OR_LOWER)
#endif
        me_context_ptr->update_hme_search_center_flag = 0;

    return NULL;
};

#if ME_HME_PRUNING_CLEANUP
void tf_set_me_hme_ref_prune_ctrls(MeContext* context_ptr, uint8_t prune_level) {
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

void tf_set_me_sr_adjustment_ctrls(MeContext* context_ptr, uint8_t sr_adjustment_level) {
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
#endif
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
#if !REFACTOR_ME_HME ||  !MAR12_ADOPTIONS
    uint8_t  hme_me_level = scs_ptr->use_output_stat_file ?
        pcs_ptr->snd_pass_enc_mode : pcs_ptr->enc_mode;
#endif
    // Set ME/HME search regions
    tf_set_me_hme_params_oq(
        context_ptr->me_context_ptr, pcs_ptr, scs_ptr, input_resolution);
#if !REFACTOR_ME_HME
    context_ptr->me_context_ptr->max_me_search_width =
        max_metf_search_width[sc_content_detected][input_resolution][hme_me_level];
    context_ptr->me_context_ptr->max_me_search_height =
        max_metf_search_height[sc_content_detected][input_resolution][hme_me_level];
#endif
#if !REMOVE_ME_SUBPEL_CODE
    if (sc_content_detected)
#if MAR11_ADOPTIONS
        // fractional_search_method is irrelevant if subpel is OFF
        context_ptr->me_context_ptr->fractional_search_method = FULL_SAD_SEARCH;
#else
        if (enc_mode <= ENC_M1)
            context_ptr->me_context_ptr->fractional_search_method = FULL_SAD_SEARCH;
        else
            context_ptr->me_context_ptr->fractional_search_method = SUB_SAD_SEARCH;
#endif
    else
#if MAR2_M8_ADOPTIONS
        context_ptr->me_context_ptr->fractional_search_method = SSD_SEARCH;
#else
    if (enc_mode <= ENC_M6)
        context_ptr->me_context_ptr->fractional_search_method = SSD_SEARCH;
    else
        context_ptr->me_context_ptr->fractional_search_method = FULL_SAD_SEARCH;
#endif
#endif
    // Set HME flags
    context_ptr->me_context_ptr->enable_hme_flag        = pcs_ptr->tf_enable_hme_flag;
    context_ptr->me_context_ptr->enable_hme_level0_flag = pcs_ptr->tf_enable_hme_level0_flag;
    context_ptr->me_context_ptr->enable_hme_level1_flag = pcs_ptr->tf_enable_hme_level1_flag;
    context_ptr->me_context_ptr->enable_hme_level2_flag = pcs_ptr->tf_enable_hme_level2_flag;
#if !REMOVE_ME_SUBPEL_CODE
    if (scs_ptr->static_config.enable_subpel == DEFAULT)
        // Set the default settings of subpel
        if (sc_content_detected)
#if MAR11_ADOPTIONS
            context_ptr->me_context_ptr->use_subpel_flag = 0;
#else
                context_ptr->me_context_ptr->use_subpel_flag = 0;
#endif
        else
            context_ptr->me_context_ptr->use_subpel_flag = 1;
    else
        context_ptr->me_context_ptr->use_subpel_flag = scs_ptr->static_config.enable_subpel;
#if M1_C2_ADOPTIONS
    if (enc_mode <= ENC_M0) {
#else
#if APR23_ADOPTIONS
    if (enc_mode <= ENC_M2) {
#else
#if MAR10_ADOPTIONS
    if (enc_mode <= ENC_M1) {
#else
    if (enc_mode <= ENC_M0) {
#endif
#endif
#endif
        context_ptr->me_context_ptr->half_pel_mode =
            (sc_content_detected) ? REFINEMENT_HP_MODE : EX_HP_MODE;
    }
#if !MAR11_ADOPTIONS
#if MAR3_M2_ADOPTIONS
    else if (enc_mode <= ENC_M2) {
#else
    else if (enc_mode <= ENC_M1) {
#endif
        context_ptr->me_context_ptr->half_pel_mode =
            (sc_content_detected) ? REFINEMENT_HP_MODE : SWITCHABLE_HP_MODE;
    }
#endif
    else {
        context_ptr->me_context_ptr->half_pel_mode    = REFINEMENT_HP_MODE;
    }
    context_ptr->me_context_ptr->h_pel_search_wind =   H_PEL_SEARCH_WIND_3;
    // Set fractional search model
    // 0: search all blocks
    // 1: selective based on Full-Search SAD & MV.
    // 2: off
    if (context_ptr->me_context_ptr->use_subpel_flag == 1) {
#if MAR17_ADOPTIONS
        context_ptr->me_context_ptr->fractional_search_model = 0;
#else
        if (enc_mode <= ENC_M6)
            context_ptr->me_context_ptr->fractional_search_model = 0;
        else
            context_ptr->me_context_ptr->fractional_search_model = 1;
#endif
    } else
        context_ptr->me_context_ptr->fractional_search_model = 2;
#endif
    // HME Search Method
#if !UNIFY_SC_NSC
    if (sc_content_detected)
#if MAR17_ADOPTIONS
        context_ptr->me_context_ptr->hme_search_method = FULL_SAD_SEARCH;
#else
        if (enc_mode <= ENC_M6)
            context_ptr->me_context_ptr->hme_search_method = FULL_SAD_SEARCH;
        else
            context_ptr->me_context_ptr->hme_search_method = SUB_SAD_SEARCH;
#endif
    else
#endif
        context_ptr->me_context_ptr->hme_search_method = FULL_SAD_SEARCH;
    // ME Search Method
#if !UNIFY_SC_NSC
    if (sc_content_detected)
#if MAR10_ADOPTIONS
        if (enc_mode <= ENC_M8)
#else
        if (enc_mode <= ENC_M3)
#endif
            context_ptr->me_context_ptr->me_search_method = FULL_SAD_SEARCH;
        else
            context_ptr->me_context_ptr->me_search_method = SUB_SAD_SEARCH;
    else
#endif
#if MAR10_ADOPTIONS
        context_ptr->me_context_ptr->me_search_method = FULL_SAD_SEARCH;
#else
        context_ptr->me_context_ptr->me_search_method =
            (enc_mode <= ENC_M4) ? FULL_SAD_SEARCH : SUB_SAD_SEARCH;
#endif
#if !SHUT_ME_NSQ_SEARCH
    // Me nsq search levels.
    // 0: feature off -> perform nsq_search.
    // 1: perform me nsq_search for the best refrenece picture.
    // 2: perform me nsq_search for the nearest refrenece pictures.
    // 3: me nsq_search off.
    context_ptr->me_context_ptr->inherit_rec_mv_from_sq_block = 0;
#endif
#if ME_HME_PRUNING_CLEANUP

    // Set hme/me based reference pruning level (0-4)
    // Ref pruning is disallowed for TF in motion_estimate_sb()
    tf_set_me_hme_ref_prune_ctrls(context_ptr->me_context_ptr, 0);

    // Set hme-based me sr adjustment level
    // ME SR adjustment is disallowed for TF in motion_estimate_sb()
    tf_set_me_sr_adjustment_ctrls(context_ptr->me_context_ptr, 0);
#else
#if ADD_ME_SIGNAL_FOR_PRUNING_TH
#if MAR20_ADOPTIONS
    if (pcs_ptr->sc_content_detected)
        context_ptr->me_context_ptr->prune_ref_if_hme_sad_dev_bigger_than_th = 200;
    else if (enc_mode <= ENC_M1)
        context_ptr->me_context_ptr->prune_ref_if_hme_sad_dev_bigger_than_th = 80;
#if PRESETS_SHIFT
    else if (enc_mode <= ENC_M2)
#else
    else if (enc_mode <= ENC_M3)
#endif
        context_ptr->me_context_ptr->prune_ref_if_hme_sad_dev_bigger_than_th = 50;
#else
#if MAR19_ADOPTIONS
    if (enc_mode <= ENC_M4)
#else
    if (enc_mode <= ENC_M7)
#endif
        context_ptr->me_context_ptr->prune_ref_if_hme_sad_dev_bigger_than_th = 80;
#endif
    else
        context_ptr->me_context_ptr->prune_ref_if_hme_sad_dev_bigger_than_th = 30;

#if MAR19_ADOPTIONS
    // Turn this threshold off
    context_ptr->me_context_ptr->prune_ref_if_me_sad_dev_bigger_than_th = (uint16_t)~0;
#else
    if (enc_mode <= ENC_M7)
        context_ptr->me_context_ptr->prune_ref_if_me_sad_dev_bigger_than_th = 30;
    else
        context_ptr->me_context_ptr->prune_ref_if_me_sad_dev_bigger_than_th = 15;
#endif
#endif
#endif
#if !ADD_MAX_HME_SIGNAL
#if ADD_HME_MIN_MAX_MULTIPLIER_SIGNAL
    context_ptr->me_context_ptr->max_hme_sr_area_multipler = 3;
#endif
#endif
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

    EbObjectWrapper *       in_results_wrapper_ptr;

    EbObjectWrapper *        out_results_wrapper_ptr;

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

        PictureDecisionResults *in_results_ptr = (PictureDecisionResults *)
                                                     in_results_wrapper_ptr->object_ptr;
        PictureParentControlSet *pcs_ptr = (PictureParentControlSet *)
                                               in_results_ptr->pcs_wrapper_ptr->object_ptr;
        SequenceControlSet *scs_ptr = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;

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
        EbPictureBufferDesc *input_padded_picture_ptr = pa_ref_obj_->input_padded_picture_ptr;

        EbPictureBufferDesc *input_picture_ptr = pcs_ptr->enhanced_picture_ptr;

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

            // Global motion estimation
            // Compute only for the first fragment.
            // TODO: create an other kernel ?
            if (pcs_ptr->gm_level == GM_FULL || pcs_ptr->gm_level == GM_DOWN) {
                if (context_ptr->me_context_ptr->compute_global_motion &&
                    in_results_ptr->segment_index == 0)
                    global_motion_estimation(
                        pcs_ptr, context_ptr->me_context_ptr, input_picture_ptr);
            }

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

                // Use scaled source references if resolution of the reference is different that of the input
                use_scaled_source_refs_if_needed(pcs_ptr,
                                                 input_picture_ptr,
                                                 pa_ref_obj_,
                                                 &input_padded_picture_ptr,
                                                 &quarter_picture_ptr,
                                                 &sixteenth_picture_ptr);

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
                        uint32_t buffer_index = (input_picture_ptr->origin_y + sb_origin_y) *
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

#if ENABLE_HME_AT_INC_SB
                            for (sb_row = 0; sb_row < (BLOCK_SIZE_64 >> 1); sb_row++) {
#else
                            for (sb_row = 0; sb_row < (sb_height >> 1); sb_row++) {
#endif
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
#if ENABLE_HME_AT_INC_SB
                                    for (sb_row = 0; sb_row < (BLOCK_SIZE_64 >> 2); sb_row++) {
#else
                                    for (sb_row = 0; sb_row < (sb_height >> 2); sb_row += 1) {
#endif
                                        eb_memcpy(local_ptr,
                                                  frame_ptr,
                                                  (sb_width >> 2) * sizeof(uint8_t));
                                        local_ptr += 16;
                                        frame_ptr += sixteenth_picture_ptr->stride_y;
                                    }
                                } else {
#if ENABLE_HME_AT_INC_SB
#if FIX_HME_LOAD
                                    for (sb_row = 0; sb_row < (BLOCK_SIZE_64 >> 2); sb_row += 2) {
#else
                                    for (sb_row = 0; sb_row < (BLOCK_SIZE_64 >> 2); sb_row++) {
#endif
#else
                                    for (sb_row = 0; sb_row < (sb_height >> 2); sb_row += 2) {
#endif
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
#if !REMOVE_UNUSED_CODE_PH2
            if (pcs_ptr->intra_pred_mode > 4)
            // *** OPEN LOOP INTRA CANDIDATE SEARCH CODE ***
            {
                // SB Loop
                for (y_sb_index = y_sb_start_index; y_sb_index < y_sb_end_index; ++y_sb_index) {
                    for (x_sb_index = x_sb_start_index; x_sb_index < x_sb_end_index; ++x_sb_index) {
                        sb_index = (uint16_t)(x_sb_index + y_sb_index * pic_width_in_sb);

                        open_loop_intra_search_sb(
                            pcs_ptr, sb_index, context_ptr, input_picture_ptr);
                    }
                }
            }
#endif
#if TPL_LA
            if (scs_ptr->static_config.look_ahead_distance != 0 && scs_ptr->static_config.enable_tpl_la) {
                for (y_sb_index = y_sb_start_index; y_sb_index < y_sb_end_index; ++y_sb_index) {
                    for (x_sb_index = x_sb_start_index; x_sb_index < x_sb_end_index; ++x_sb_index) {
                        sb_origin_x = x_sb_index * scs_ptr->sb_sz;
                        sb_origin_y = y_sb_index * scs_ptr->sb_sz;

                        sb_index = (uint16_t)(x_sb_index + y_sb_index * pic_width_in_sb);
                        open_loop_intra_search_mb(pcs_ptr, sb_index, input_picture_ptr);
                    }
                }
            }
#endif
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

            MotionEstimationResults *out_results_ptr = (MotionEstimationResults *)
                                                           out_results_wrapper_ptr->object_ptr;
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
