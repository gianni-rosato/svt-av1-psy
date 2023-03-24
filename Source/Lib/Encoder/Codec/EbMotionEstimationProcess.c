// clang-format off
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
EbErrorType check_00_center(PictureParentControlSet *pcs, EbPictureBufferDesc *ref_pic_ptr,
                            MeContext *me_ctx, uint32_t sb_origin_x, uint32_t sb_origin_y,
                            uint32_t sb_width, uint32_t sb_height, int16_t *x_search_center,
                            int16_t *y_search_center, uint32_t zz_sad);

/************************************************
 * Set ME/HME Params from Config
 ************************************************/
/************************************************
 * Set HME Search area parameters
 ************************************************/
void set_hme_search_params(PictureParentControlSet *pcs, MeContext *me_ctx,
                           EbInputResolution input_resolution) {
    // Set number of HME level 0 search regions to use
    me_ctx->num_hme_sa_w = 2;
    me_ctx->num_hme_sa_h = 2;

    // Set HME level 0 min and max search areas
    if (pcs->enc_mode <= ENC_MRS) {
        me_ctx->hme_l0_sa.sa_min = (SearchArea){ 240, 240 };
        me_ctx->hme_l0_sa.sa_max = (SearchArea){ 480, 480 };
    } else if (pcs->enc_mode <= ENC_M0) {
        if (input_resolution < INPUT_SIZE_4K_RANGE) {
            me_ctx->hme_l0_sa.sa_min = (SearchArea){32, 32};
            me_ctx->hme_l0_sa.sa_max = (SearchArea){192, 192};
        } else {
            me_ctx->hme_l0_sa.sa_min = (SearchArea){240, 240};
            me_ctx->hme_l0_sa.sa_max = (SearchArea){480, 480};
        }
    }
    else if (pcs->enc_mode <= ENC_M4) {
        me_ctx->hme_l0_sa.sa_min = (SearchArea){32, 32};
        me_ctx->hme_l0_sa.sa_max = (SearchArea){192, 192};
    }
    else if (pcs->enc_mode <= ENC_M6) {
        if (input_resolution <= INPUT_SIZE_1080p_RANGE) {
            if (pcs->sc_class1) {
                me_ctx->hme_l0_sa.sa_min = (SearchArea){32, 32};
                me_ctx->hme_l0_sa.sa_max = (SearchArea){192, 192};
            } else {
                me_ctx->hme_l0_sa.sa_min = (SearchArea) { 16, 16 };
                me_ctx->hme_l0_sa.sa_max = (SearchArea) { 192, 192 };
            }
        } else {
            me_ctx->hme_l0_sa.sa_min = (SearchArea){32, 32};
            me_ctx->hme_l0_sa.sa_max = (SearchArea){192, 192};
        }
    }

    else if (pcs->enc_mode <= ENC_M11) {
        if (pcs->sc_class1) {
            me_ctx->hme_l0_sa.sa_min = (SearchArea){32, 32};
            me_ctx->hme_l0_sa.sa_max = (SearchArea){192, 192};
        } else {
            me_ctx->hme_l0_sa.sa_min = (SearchArea) { 16, 16 };
            me_ctx->hme_l0_sa.sa_max = (SearchArea) { 192, 192 };
        }
    }
    else {
        if (pcs->sc_class1) {
            me_ctx->hme_l0_sa.sa_min = (SearchArea){ 32, 32 };
            me_ctx->hme_l0_sa.sa_max = (SearchArea){ 192, 192 };
        }
        else {
            if (input_resolution < INPUT_SIZE_4K_RANGE) {
                me_ctx->hme_l0_sa.sa_min = (SearchArea){ 8, 8 };
                me_ctx->hme_l0_sa.sa_max = (SearchArea){ 96, 96 };
            }
            else {
                me_ctx->hme_l0_sa.sa_min = (SearchArea){ 16, 16 };
                me_ctx->hme_l0_sa.sa_max = (SearchArea){ 96, 96 };
            }
        }
    }
    // Set the HME Level 1 and Level 2 refinement areas
    if (pcs->enc_mode <= ENC_M1) {
        me_ctx->hme_l1_sa = (SearchArea){16, 16};
        me_ctx->hme_l2_sa = (SearchArea){16, 16};
    } else {
        me_ctx->hme_l1_sa = (SearchArea){8, 3};
        me_ctx->hme_l2_sa = (SearchArea){8, 3};
    }
}

/************************************************
 * Set ME Search area parameters
 ************************************************/
void set_me_search_params(SequenceControlSet *scs, PictureParentControlSet *pcs,
                          MeContext *me_ctx, EbInputResolution input_resolution) {
#if FIX_LAYER_SIGNAL
    uint32_t hierarchical_levels = pcs->hierarchical_levels;
#else
    uint32_t hierarchical_levels = scs->static_config.hierarchical_levels;
#endif
    // Set the min and max ME search area
    if (pcs->sc_class1) {
        if (pcs->enc_mode <= ENC_M2) {
            me_ctx->me_sa.sa_min = (SearchArea){175, 175};
            me_ctx->me_sa.sa_max = (SearchArea){750, 750};
        } else if (pcs->enc_mode <= ENC_M4) {
            me_ctx->me_sa.sa_min = (SearchArea){125, 125};
            me_ctx->me_sa.sa_max = (SearchArea){500, 500};
        } else if (pcs->enc_mode <= ENC_M6) {
            me_ctx->me_sa.sa_min = (SearchArea){75, 75};
            me_ctx->me_sa.sa_max = (SearchArea){350, 350};
        } else if (pcs->enc_mode <= ENC_M7) {
            me_ctx->me_sa.sa_min = (SearchArea){48, 48};
            me_ctx->me_sa.sa_max = (SearchArea){224, 224};
        } else if (pcs->enc_mode <= ENC_M8) {
            me_ctx->me_sa.sa_min = (SearchArea){32, 32};
            me_ctx->me_sa.sa_max = (SearchArea){164, 164};
        } else if (pcs->enc_mode <= ENC_M10) {
            me_ctx->me_sa.sa_min = (SearchArea){32, 32};
            me_ctx->me_sa.sa_max = (SearchArea){96, 96};
        } else if (pcs->enc_mode <= ENC_M12) {
            me_ctx->me_sa.sa_min = (SearchArea) { 16, 16 };
            me_ctx->me_sa.sa_max = (SearchArea) { 96, 96 };
        }
        else {
            me_ctx->me_sa.sa_min = (SearchArea) { 8, 8 };
            me_ctx->me_sa.sa_max = (SearchArea) { 32, 32 };
        }
    } else if (pcs->enc_mode <= ENC_M0) {
        me_ctx->me_sa.sa_min = (SearchArea){64, 64};
        me_ctx->me_sa.sa_max = (SearchArea){256, 256};
    } else if (pcs->enc_mode <= ENC_M2) {
        me_ctx->me_sa.sa_min = (SearchArea){64, 64};
        me_ctx->me_sa.sa_max = (SearchArea){128, 128};
    }
    else if (pcs->enc_mode <= ENC_M3) {
        me_ctx->me_sa.sa_min = (SearchArea){16, 16};
        me_ctx->me_sa.sa_max = (SearchArea){128, 128};
    }
    else if (pcs->enc_mode <= ENC_M4) {
        me_ctx->me_sa.sa_min = (SearchArea){16, 16};
        me_ctx->me_sa.sa_max = (SearchArea){64, 64};
    }
    else if (pcs->enc_mode <= ENC_M5) {
        me_ctx->me_sa.sa_min = (SearchArea) { 16, 16 };
        me_ctx->me_sa.sa_max = (SearchArea) { 64, 32 };
    }
    else if (pcs->enc_mode <= ENC_M7) {
        if (input_resolution < INPUT_SIZE_1080p_RANGE) {
            me_ctx->me_sa.sa_min = (SearchArea) { 16, 16 };
            me_ctx->me_sa.sa_max = (SearchArea) { 32, 16 };
        }
        else {
            me_ctx->me_sa.sa_min = (SearchArea) { 16, 6 };
            me_ctx->me_sa.sa_max = (SearchArea) { 16, 9 };
        }
    }
    else if (pcs->enc_mode <= ENC_M8) {
        if (hierarchical_levels <= 3) {
            if (input_resolution < INPUT_SIZE_4K_RANGE) {
                me_ctx->me_sa.sa_min = (SearchArea) { 8, 5 };
                me_ctx->me_sa.sa_max = (SearchArea) { 16, 9 };
            }
            else {
                me_ctx->me_sa.sa_min = (SearchArea) { 8, 1 };
                me_ctx->me_sa.sa_max = (SearchArea) { 8, 1 };
            }
        }
        else if (input_resolution < INPUT_SIZE_1080p_RANGE) {
            me_ctx->me_sa.sa_min = (SearchArea) { 16, 16 };
            me_ctx->me_sa.sa_max = (SearchArea) { 32, 16 };
        }
        else {
            me_ctx->me_sa.sa_min = (SearchArea) { 16, 6 };
            me_ctx->me_sa.sa_max = (SearchArea) { 16, 9 };
        }
    }
    else if (pcs->enc_mode <= ENC_M9) {
        if (input_resolution < INPUT_SIZE_1080p_RANGE) {
            me_ctx->me_sa.sa_min = (SearchArea) { 16, 16 };
            me_ctx->me_sa.sa_max = (SearchArea) { 32, 16 };
        }
        else {
            me_ctx->me_sa.sa_min = (SearchArea) { 16, 6 };
            me_ctx->me_sa.sa_max = (SearchArea) { 16, 9 };
        }
    }
    else if (pcs->enc_mode <= ENC_M10) {
        if (input_resolution < INPUT_SIZE_4K_RANGE) {
            me_ctx->me_sa.sa_min = (SearchArea) { 8, 5 };
            me_ctx->me_sa.sa_max = (SearchArea) { 16, 9 };
        }
        else {
            me_ctx->me_sa.sa_min = (SearchArea) { 8, 1 };
            me_ctx->me_sa.sa_max = (SearchArea) { 8, 1 };
        }
    }
    else if (pcs->enc_mode <= ENC_M12) {
        if (input_resolution < INPUT_SIZE_720p_RANGE) {
            me_ctx->me_sa.sa_min = (SearchArea){8, 3};
            me_ctx->me_sa.sa_max = (SearchArea){16, 9};
        } else if (input_resolution < INPUT_SIZE_1080p_RANGE) {
            me_ctx->me_sa.sa_min = (SearchArea){8, 1};
            me_ctx->me_sa.sa_max = (SearchArea){16, 7};
        } else if (input_resolution < INPUT_SIZE_4K_RANGE) {
            me_ctx->me_sa.sa_min = (SearchArea){8, 1};
            me_ctx->me_sa.sa_max = (SearchArea){8, 7};
        } else {
            me_ctx->me_sa.sa_min = (SearchArea){8, 1};
            me_ctx->me_sa.sa_max = (SearchArea){8, 1};
        }
    } else {
        me_ctx->me_sa.sa_min = (SearchArea){8, 3};
        me_ctx->me_sa.sa_max = (SearchArea){8, 3};
    }
#if OPT_LD_M9
    const bool               rtc_tune = (scs->static_config.pred_structure == SVT_AV1_PRED_LOW_DELAY_B) ? true : false;
    if (rtc_tune && pcs->enc_mode > ENC_M8) {
        if (input_resolution < INPUT_SIZE_720p_RANGE) {
            me_ctx->me_sa.sa_min = (SearchArea) { 8, 3 };
            me_ctx->me_sa.sa_max = (SearchArea) { 16, 9 };
        }
        else if (input_resolution < INPUT_SIZE_1080p_RANGE) {
            me_ctx->me_sa.sa_min = (SearchArea) { 8, 1 };
            me_ctx->me_sa.sa_max = (SearchArea) { 16, 7 };
        }
        else if (input_resolution < INPUT_SIZE_4K_RANGE) {
            me_ctx->me_sa.sa_min = (SearchArea) { 8, 1 };
            me_ctx->me_sa.sa_max = (SearchArea) { 8, 7 };
        }
        else {
            me_ctx->me_sa.sa_min = (SearchArea) { 8, 1 };
            me_ctx->me_sa.sa_max = (SearchArea) { 8, 1 };
        }
    }
#endif
    // Scale up the MIN ME area if low frame rate
    bool low_frame_rate_flag = (scs->frame_rate >> 16);
    if (low_frame_rate_flag) {
        me_ctx->me_sa.sa_min.width  = (me_ctx->me_sa.sa_min.width * 3) >> 1;
        me_ctx->me_sa.sa_min.height = (me_ctx->me_sa.sa_min.height * 3) >> 1;
    }
}
void svt_aom_set_me_hme_ref_prune_ctrls(MeContext *me_ctx, uint8_t prune_level) {
    MeHmeRefPruneCtrls *me_hme_prune_ctrls = &me_ctx->me_hme_prune_ctrls;

    switch (prune_level) {
    case 0:
        me_hme_prune_ctrls->enable_me_hme_ref_pruning               = 0;
        me_hme_prune_ctrls->prune_ref_if_hme_sad_dev_bigger_than_th = (uint16_t)~0;
        me_hme_prune_ctrls->prune_ref_if_me_sad_dev_bigger_than_th  = (uint16_t)~0;
        break;
    case 1:
        me_hme_prune_ctrls->enable_me_hme_ref_pruning               = 1;
        me_hme_prune_ctrls->prune_ref_if_hme_sad_dev_bigger_than_th = 160;
        me_hme_prune_ctrls->prune_ref_if_me_sad_dev_bigger_than_th  = (uint16_t)~0;
        me_hme_prune_ctrls->protect_closest_refs                    = 1;
        break;
    case 2:
        me_hme_prune_ctrls->enable_me_hme_ref_pruning               = 1;
        me_hme_prune_ctrls->prune_ref_if_hme_sad_dev_bigger_than_th = 80;
        me_hme_prune_ctrls->prune_ref_if_me_sad_dev_bigger_than_th  = 60;
        me_hme_prune_ctrls->protect_closest_refs                    = 1;
        break;
    case 3:
        me_hme_prune_ctrls->enable_me_hme_ref_pruning               = 1;
        me_hme_prune_ctrls->prune_ref_if_hme_sad_dev_bigger_than_th = 50;
        me_hme_prune_ctrls->prune_ref_if_me_sad_dev_bigger_than_th  = 60;
        me_hme_prune_ctrls->protect_closest_refs                    = 1;
        break;
    case 4:
        me_hme_prune_ctrls->enable_me_hme_ref_pruning               = 1;
        me_hme_prune_ctrls->prune_ref_if_hme_sad_dev_bigger_than_th = 30;
        me_hme_prune_ctrls->prune_ref_if_me_sad_dev_bigger_than_th  = 60;
        me_hme_prune_ctrls->protect_closest_refs                    = 1;
        break;
    case 5:
        me_hme_prune_ctrls->enable_me_hme_ref_pruning               = 1;
        me_hme_prune_ctrls->prune_ref_if_hme_sad_dev_bigger_than_th = 5;
        me_hme_prune_ctrls->prune_ref_if_me_sad_dev_bigger_than_th  = 60;
        me_hme_prune_ctrls->protect_closest_refs                    = 1;
        break;
    case 6:
        me_hme_prune_ctrls->enable_me_hme_ref_pruning               = 1;
        me_hme_prune_ctrls->prune_ref_if_hme_sad_dev_bigger_than_th = 200;
        me_hme_prune_ctrls->prune_ref_if_me_sad_dev_bigger_than_th  = 60;
        me_hme_prune_ctrls->protect_closest_refs                    = 0;
        break;
    default: assert(0); break;
    }
}

void svt_aom_set_me_sr_adjustment_ctrls(MeContext *me_ctx, uint8_t sr_adjustment_level) {
    MeSrCtrls *me_sr_adjustment_ctrls = &me_ctx->me_sr_adjustment_ctrls;

    switch (sr_adjustment_level) {
    case 0: me_sr_adjustment_ctrls->enable_me_sr_adjustment = 0; break;
    case 1:
        me_sr_adjustment_ctrls->enable_me_sr_adjustment              = 1;
        me_sr_adjustment_ctrls->reduce_me_sr_based_on_mv_length_th   = 4;
        me_sr_adjustment_ctrls->stationary_hme_sad_abs_th            = 12000;
        me_sr_adjustment_ctrls->stationary_me_sr_divisor             = 8;
        me_sr_adjustment_ctrls->reduce_me_sr_based_on_hme_sad_abs_th = 6000;
        me_sr_adjustment_ctrls->me_sr_divisor_for_low_hme_sad        = 8;
        me_sr_adjustment_ctrls->distance_based_hme_resizing          = 0;
        break;
    case 2:
        me_sr_adjustment_ctrls->enable_me_sr_adjustment              = 1;
        me_sr_adjustment_ctrls->reduce_me_sr_based_on_mv_length_th   = 4;
        me_sr_adjustment_ctrls->stationary_hme_sad_abs_th            = 12000;
        me_sr_adjustment_ctrls->stationary_me_sr_divisor             = 8;
        me_sr_adjustment_ctrls->reduce_me_sr_based_on_hme_sad_abs_th = 6000;
        me_sr_adjustment_ctrls->me_sr_divisor_for_low_hme_sad        = 8;
        me_sr_adjustment_ctrls->distance_based_hme_resizing          = 1;
        break;
    case 3:
        me_sr_adjustment_ctrls->enable_me_sr_adjustment              = 1;
        me_sr_adjustment_ctrls->reduce_me_sr_based_on_mv_length_th   = 4;
        me_sr_adjustment_ctrls->stationary_hme_sad_abs_th            = 12000;
        me_sr_adjustment_ctrls->stationary_me_sr_divisor             = 8;
        me_sr_adjustment_ctrls->reduce_me_sr_based_on_hme_sad_abs_th = 12000;
        me_sr_adjustment_ctrls->me_sr_divisor_for_low_hme_sad        = 8;
        me_sr_adjustment_ctrls->distance_based_hme_resizing          = 1;
        break;
    case 4:
        me_sr_adjustment_ctrls->enable_me_sr_adjustment = 2;
        me_sr_adjustment_ctrls->reduce_me_sr_based_on_mv_length_th = 16;
        me_sr_adjustment_ctrls->stationary_hme_sad_abs_th = 20000;
        me_sr_adjustment_ctrls->stationary_me_sr_divisor = 8;
        me_sr_adjustment_ctrls->reduce_me_sr_based_on_hme_sad_abs_th = 20000;
        me_sr_adjustment_ctrls->me_sr_divisor_for_low_hme_sad = 8;
        me_sr_adjustment_ctrls->distance_based_hme_resizing   = 1;
        break;

    case 5:

        me_sr_adjustment_ctrls->enable_me_sr_adjustment = 2;
        me_sr_adjustment_ctrls->reduce_me_sr_based_on_mv_length_th = 20;
        me_sr_adjustment_ctrls->stationary_hme_sad_abs_th = 24000;
        me_sr_adjustment_ctrls->stationary_me_sr_divisor = 8;
        me_sr_adjustment_ctrls->reduce_me_sr_based_on_hme_sad_abs_th = 24000;
        me_sr_adjustment_ctrls->me_sr_divisor_for_low_hme_sad = 8;
        me_sr_adjustment_ctrls->distance_based_hme_resizing   = 1;

        break;

    default: assert(0); break;
    }
    if (me_ctx->enable_hme_level2_flag == 0) {
        if (me_ctx->enable_hme_level1_flag == 1) {
            me_sr_adjustment_ctrls->stationary_hme_sad_abs_th =
                me_sr_adjustment_ctrls->stationary_hme_sad_abs_th / 4;
            me_sr_adjustment_ctrls->reduce_me_sr_based_on_hme_sad_abs_th =
                me_sr_adjustment_ctrls->reduce_me_sr_based_on_hme_sad_abs_th / 4;
        } else {
            me_sr_adjustment_ctrls->stationary_hme_sad_abs_th =
                me_sr_adjustment_ctrls->stationary_hme_sad_abs_th / 16;
            me_sr_adjustment_ctrls->reduce_me_sr_based_on_hme_sad_abs_th =
                me_sr_adjustment_ctrls->reduce_me_sr_based_on_hme_sad_abs_th / 16;
        }
    }
}
/*set the skip_frame flag for ipp*/
void set_skip_frame_in_ipp(PictureParentControlSet * pcs,MeContext *me_ctx) {
    me_ctx->skip_frame = 0;
    int8_t max_heirachical_level = (int8_t)pcs->scs->max_heirachical_level;
    int8_t enc_mode = pcs->scs->static_config.enc_mode;
    int8_t pass = (int8_t)pcs->scs->static_config.pass;
    uint8_t skip_frame_first_pass = pcs->scs->ipp_pass_ctrls.skip_frame_first_pass;
    uint64_t  picture_number = pcs->picture_number;
    if (!skip_frame_first_pass) {
        if (enc_mode < ENC_M8)
            me_ctx->skip_frame = 0;
        else if (pass == ENC_SINGLE_PASS)
        {
            // SINGLE PASS
            if (max_heirachical_level <= 2)
                me_ctx->skip_frame = ((enc_mode >= ENC_M10) && (picture_number > 1 && picture_number % 2 > 0)) ? 1 : 0;
            else if ((enc_mode < ENC_M10 && (picture_number > 3 && picture_number % 4 > 0)) ||
                (enc_mode >= ENC_M10 && (picture_number > 5 && picture_number % 6 > 0))) {
                me_ctx->skip_frame = 1;
            }
        }
    }
    else {
        if ((skip_frame_first_pass == 1) &&
            (picture_number % 8 > 0))
            me_ctx->skip_frame = 1;
        else if ((skip_frame_first_pass == 2) &&
            (picture_number > 7 && picture_number % 8 > 0))
            me_ctx->skip_frame = 1;
        else
            if (picture_number > 3 && picture_number % 4 > 0)
                me_ctx->skip_frame = 1;
            else
                me_ctx->skip_frame = 0;
    }
}
/*configure PreHme control*/
void svt_aom_set_prehme_ctrls(MeContext *me_ctx, uint8_t level) {
    PreHmeCtrls *ctrl = &me_ctx->prehme_ctrl;

    switch (level) {
    case 0: ctrl->enable = 0; break;
    case 1:
        ctrl->enable = 1;
        // vertical shape search region
        ctrl->prehme_sa_cfg[0].sa_min = (SearchArea){8, 100};
        ctrl->prehme_sa_cfg[0].sa_max = (SearchArea){8, 400};
        // horizontal shape search region
        ctrl->prehme_sa_cfg[1].sa_min = (SearchArea){96, 3};
        ctrl->prehme_sa_cfg[1].sa_max = (SearchArea){384, 3};
        ctrl->skip_search_line = 0;
        ctrl->l1_early_exit = 0;
        break;
    case 2: ctrl->enable = 1;
        // vertical shape search region
        ctrl->prehme_sa_cfg[0].sa_min = (SearchArea){8, 100};
        ctrl->prehme_sa_cfg[0].sa_max = (SearchArea){8, 350};
        // horizontal shape search region
        ctrl->prehme_sa_cfg[1].sa_min = (SearchArea){32, 7};
        ctrl->prehme_sa_cfg[1].sa_max = (SearchArea){200, 7};
        ctrl->skip_search_line = 1;
        ctrl->l1_early_exit = 0;
        break;
    case 3:
        ctrl->enable = 1;
        // vertical shape search region
        ctrl->prehme_sa_cfg[0].sa_min = (SearchArea){8, 100};
        ctrl->prehme_sa_cfg[0].sa_max = (SearchArea){8, 350};
        // horizontal shape search region
        ctrl->prehme_sa_cfg[1].sa_min = (SearchArea){32, 7};
        ctrl->prehme_sa_cfg[1].sa_max = (SearchArea){128, 7};
        ctrl->skip_search_line        = 1;
        ctrl->l1_early_exit = 1;
        break;
    default: assert(0); break;
    }
}
/******************************************************
* Derive ME Settings for OQ
  Input   : encoder mode and tune
  Output  : ME Kernel signal(s)
******************************************************/
EbErrorType svt_aom_signal_derivation_me_kernel_oq(SequenceControlSet *       scs,
                                           PictureParentControlSet *  pcs,
                                           MotionEstimationContext_t *me_context_ptr) {
    EbErrorType return_error = EB_ErrorNone;

    EncMode         enc_mode         = pcs->enc_mode;
    EbInputResolution input_resolution = scs->input_resolution;
#if FIX_LAYER_SIGNAL
    const uint32_t hierarchical_levels = pcs->hierarchical_levels;
#else
    const uint32_t hierarchical_levels = scs->static_config.hierarchical_levels;
#endif
    // Set ME search area
    set_me_search_params(scs, pcs, me_context_ptr->me_ctx, input_resolution);

    // Set HME search area
    set_hme_search_params(pcs, me_context_ptr->me_ctx, input_resolution);

    // Set HME flags
    me_context_ptr->me_ctx->enable_hme_flag        = pcs->enable_hme_flag;
    me_context_ptr->me_ctx->enable_hme_level0_flag = pcs->enable_hme_level0_flag;
    me_context_ptr->me_ctx->enable_hme_level1_flag = pcs->enable_hme_level1_flag;
    me_context_ptr->me_ctx->enable_hme_level2_flag = pcs->enable_hme_level2_flag;
    // HME Search Method
    me_context_ptr->me_ctx->hme_search_method = SUB_SAD_SEARCH;
    me_context_ptr->me_ctx->me_search_method  = SUB_SAD_SEARCH;


    if (pcs->sc_class1) {
        me_context_ptr->me_ctx->reduce_hme_l0_sr_th_min = 0;
        me_context_ptr->me_ctx->reduce_hme_l0_sr_th_max = 0;
    }
    else if (pcs->enc_mode <= ENC_M10) {
        me_context_ptr->me_ctx->reduce_hme_l0_sr_th_min = 0;
        me_context_ptr->me_ctx->reduce_hme_l0_sr_th_max = 0;
    }
    else {
        me_context_ptr->me_ctx->reduce_hme_l0_sr_th_min = 8;
        me_context_ptr->me_ctx->reduce_hme_l0_sr_th_max = 200;
    }
    // Set pre-hme level (0-2)
    uint8_t prehme_level = 0;
    if (pcs->sc_class1)
        prehme_level = 1;
    else
    {
        if (enc_mode <= ENC_M10)
            prehme_level = 1;
        else
            prehme_level = 3;
    }
    if (pcs->enable_hme_level1_flag == 0)
        prehme_level = 0;

    svt_aom_set_prehme_ctrls(me_context_ptr->me_ctx, prehme_level);

    // Set hme/me based reference pruning level (0-4)
    if (pcs->sc_class1) {
        if (enc_mode <= ENC_MRS)
            svt_aom_set_me_hme_ref_prune_ctrls(me_context_ptr->me_ctx, 0);
        else if (enc_mode <= ENC_M1)
            svt_aom_set_me_hme_ref_prune_ctrls(me_context_ptr->me_ctx, 2);
        else if (enc_mode <= ENC_M9)
            svt_aom_set_me_hme_ref_prune_ctrls(me_context_ptr->me_ctx, 5);
        else
            svt_aom_set_me_hme_ref_prune_ctrls(me_context_ptr->me_ctx, 6);
    } else {
        if (enc_mode <= ENC_MRS)
            svt_aom_set_me_hme_ref_prune_ctrls(me_context_ptr->me_ctx, 0);
        else if (enc_mode <= ENC_M0)
            svt_aom_set_me_hme_ref_prune_ctrls(me_context_ptr->me_ctx, 2);
        else if (enc_mode <= ENC_M1) {
            if (pcs->temporal_layer_index == 0)
                svt_aom_set_me_hme_ref_prune_ctrls(me_context_ptr->me_ctx, 2);
            else
                svt_aom_set_me_hme_ref_prune_ctrls(me_context_ptr->me_ctx, 4);
        }
        else if (enc_mode <= ENC_M7) {
            if (pcs->temporal_layer_index == 0)
                svt_aom_set_me_hme_ref_prune_ctrls(me_context_ptr->me_ctx, 2);
            else
                svt_aom_set_me_hme_ref_prune_ctrls(me_context_ptr->me_ctx, 5);
        }
        else
            svt_aom_set_me_hme_ref_prune_ctrls(me_context_ptr->me_ctx, 5);
    }
    // Set hme-based me sr adjustment level
    if (pcs->sc_class1)
        if (enc_mode <= ENC_M9)
            svt_aom_set_me_sr_adjustment_ctrls(me_context_ptr->me_ctx, 4);
        else
            svt_aom_set_me_sr_adjustment_ctrls(me_context_ptr->me_ctx, 5);
    else if (enc_mode <= ENC_MRS)
        svt_aom_set_me_sr_adjustment_ctrls(me_context_ptr->me_ctx, 0);
    else if (enc_mode <= ENC_M2)
        svt_aom_set_me_sr_adjustment_ctrls(me_context_ptr->me_ctx, 1);
    else
        svt_aom_set_me_sr_adjustment_ctrls(me_context_ptr->me_ctx, 3);
    if (enc_mode <= ENC_M2)
        me_context_ptr->me_ctx->prune_me_candidates_th = 0;
    else
        me_context_ptr->me_ctx->prune_me_candidates_th = 65;
    // Set signal at picture level b/c may check signal in MD
    me_context_ptr->me_ctx->use_best_unipred_cand_only =
        pcs->use_best_me_unipred_cand_only;
    if (enc_mode <= ENC_M6)
        me_context_ptr->me_ctx->me_early_exit_th = 0;
    else if (enc_mode <= ENC_M7) {
        if (hierarchical_levels <= 3)
            me_context_ptr->me_ctx->me_early_exit_th = BLOCK_SIZE_64 * BLOCK_SIZE_64 * 8;
        else
            me_context_ptr->me_ctx->me_early_exit_th = 0;
    }
    else if (pcs->enc_mode <= ENC_M10)
        me_context_ptr->me_ctx->me_early_exit_th = BLOCK_SIZE_64 * BLOCK_SIZE_64 * 8;
    else
        me_context_ptr->me_ctx->me_early_exit_th = BLOCK_SIZE_64 * BLOCK_SIZE_64 * 9;

    me_context_ptr->me_ctx->skip_frame = 0;
    return return_error;
};
void open_loop_first_pass(struct PictureParentControlSet *ppcs,
                          MotionEstimationContext_t *me_context_ptr, int32_t segment_index);

/******************************************************
* Derive ME Settings for first pass
  Input   : encoder mode and tune
  Output  : ME Kernel signal(s)
******************************************************/
EbErrorType first_pass_signal_derivation_me_kernel(SequenceControlSet *       scs,
                                                   PictureParentControlSet *  pcs,
                                                   MotionEstimationContext_t *me_context_ptr);

/************************************************
 * Set ME/HME Params for Altref Temporal Filtering
 ************************************************/
static void tf_set_me_hme_params_oq(MeContext *me_ctx, PictureParentControlSet *pcs) {

    switch (pcs->tf_ctrls.hme_me_level) {
    case 0:
        me_ctx->num_hme_sa_w = 2;
        me_ctx->num_hme_sa_h = 2;
        me_ctx->hme_l0_sa.sa_min = (SearchArea){30, 30};
        me_ctx->hme_l0_sa.sa_max = (SearchArea){60, 60};
        me_ctx->hme_l1_sa        = (SearchArea){16, 16};
        me_ctx->hme_l2_sa        = (SearchArea){16, 16};
        me_ctx->me_sa.sa_min     = (SearchArea){60, 60};
        me_ctx->me_sa.sa_max     = (SearchArea){120, 120};
        break;

    case 1:
        me_ctx->num_hme_sa_w = 2;
        me_ctx->num_hme_sa_h = 2;
        me_ctx->hme_l0_sa.sa_min = (SearchArea){16, 16};
        me_ctx->hme_l0_sa.sa_max = (SearchArea){32, 32};
        me_ctx->hme_l1_sa        = (SearchArea){16, 16};
        me_ctx->hme_l2_sa        = (SearchArea){16, 16};
        me_ctx->me_sa.sa_min     = (SearchArea){16, 16};
        me_ctx->me_sa.sa_max     = (SearchArea){32, 32};
        break;

    case 2:
        me_ctx->num_hme_sa_w = 2;
        me_ctx->num_hme_sa_h = 2;
        me_ctx->hme_l0_sa.sa_min = (SearchArea){8, 8};
        me_ctx->hme_l0_sa.sa_max = (SearchArea){16, 16};
        me_ctx->hme_l1_sa        = (SearchArea){16, 16};
        me_ctx->hme_l2_sa        = (SearchArea){16, 16};
        me_ctx->me_sa.sa_min     = (SearchArea){8, 4};
        me_ctx->me_sa.sa_max     = (SearchArea){16, 8};
        break;
#if OPT_LD_TF
    case 3:
        me_ctx->num_hme_sa_w = 2;
        me_ctx->num_hme_sa_h = 2;
        me_ctx->hme_l0_sa.sa_min = (SearchArea) { 4, 4 };
        me_ctx->hme_l0_sa.sa_max = (SearchArea) { 4, 4 };
        me_ctx->hme_l1_sa = (SearchArea) { 8, 8 };
        me_ctx->hme_l2_sa = (SearchArea) { 8, 8 };
        me_ctx->me_sa.sa_min = (SearchArea) { 8, 8 };
        me_ctx->me_sa.sa_max = (SearchArea) { 8, 8 };
        break;
#endif

    default: assert(0); break;
    }
};
/******************************************************
* Derive ME Settings for OQ for Altref Temporal Filtering
  Input   : encoder mode and tune
  Output  : ME Kernel signal(s)
******************************************************/
EbErrorType tf_signal_derivation_me_kernel_oq(PictureParentControlSet *  pcs,
                                              MotionEstimationContext_t *me_context_ptr) {
    EbErrorType return_error = EB_ErrorNone;
    // Set ME/HME search regions
    tf_set_me_hme_params_oq(me_context_ptr->me_ctx, pcs);
    // Set HME flags
    me_context_ptr->me_ctx->enable_hme_flag        = pcs->tf_enable_hme_flag;
    me_context_ptr->me_ctx->enable_hme_level0_flag = pcs->tf_enable_hme_level0_flag;
    me_context_ptr->me_ctx->enable_hme_level1_flag = pcs->tf_enable_hme_level1_flag;
    me_context_ptr->me_ctx->enable_hme_level2_flag = pcs->tf_enable_hme_level2_flag;
    // HME Search Method
    me_context_ptr->me_ctx->hme_search_method = SUB_SAD_SEARCH;
    // ME Search Method
    me_context_ptr->me_ctx->me_search_method = SUB_SAD_SEARCH;

    uint8_t prehme_level = 0;
    svt_aom_set_prehme_ctrls(me_context_ptr->me_ctx, prehme_level);

    // Set hme/me based reference pruning level (0-4)
    // Ref pruning is disallowed for TF in motion_estimate_sb()
    svt_aom_set_me_hme_ref_prune_ctrls(me_context_ptr->me_ctx, 0);

    // Set hme-based me sr adjustment level
    // ME SR adjustment is disallowed for TF in motion_estimate_sb()
    svt_aom_set_me_sr_adjustment_ctrls(me_context_ptr->me_ctx, 0);

    me_context_ptr->me_ctx->me_early_exit_th = 0;
    me_context_ptr->me_ctx->reduce_hme_l0_sr_th_min = 0;
    me_context_ptr->me_ctx->reduce_hme_l0_sr_th_max = 0;
    me_context_ptr->me_ctx->skip_frame = 0;
    return return_error;
};

static void motion_estimation_context_dctor(EbPtr p) {
    EbThreadContext *          thread_ctx = (EbThreadContext *)p;
    MotionEstimationContext_t *obj = (MotionEstimationContext_t *)thread_ctx->priv;
    EB_DELETE(obj->me_ctx);
    EB_FREE_ARRAY(obj);
}

/************************************************
 * Motion Analysis Context Constructor
 ************************************************/
EbErrorType svt_aom_motion_estimation_context_ctor(EbThreadContext *  thread_ctx,
                                           const EbEncHandle *enc_handle_ptr, int index) {
    MotionEstimationContext_t *me_context_ptr;

    EB_CALLOC_ARRAY(me_context_ptr, 1);
    thread_ctx->priv                             = me_context_ptr;
    thread_ctx->dctor                            = motion_estimation_context_dctor;
    me_context_ptr->picture_decision_results_input_fifo_ptr = svt_system_resource_get_consumer_fifo(
        enc_handle_ptr->picture_decision_results_resource_ptr, index);
    me_context_ptr->motion_estimation_results_output_fifo_ptr = svt_system_resource_get_producer_fifo(
        enc_handle_ptr->motion_estimation_results_resource_ptr, index);
    EB_NEW(me_context_ptr->me_ctx, svt_aom_me_context_ctor);
    return EB_ErrorNone;
}
/************************************************
 * Motion Analysis Kernel
 * The Motion Analysis performs  Motion Estimation
 * This process has access to the current input picture as well as
 * the input pictures, which the current picture references according
 * to the prediction structure pattern.  The Motion Analysis process is multithreaded,
 * so pictures can be processed out of order as long as all inputs are available.
 ************************************************/
void *svt_aom_motion_estimation_kernel(void *input_ptr) {
    EbThreadContext *          thread_ctx = (EbThreadContext *)input_ptr;
    MotionEstimationContext_t *me_context_ptr = (MotionEstimationContext_t *)thread_ctx->priv;
    EbObjectWrapper *          in_results_wrapper_ptr;
    EbObjectWrapper *          out_results_wrapper;
    for (;;) {
        // Get Input Full Object
        EB_GET_FULL_OBJECT(me_context_ptr->picture_decision_results_input_fifo_ptr,
                           &in_results_wrapper_ptr);
        PictureDecisionResults *in_results_ptr = (PictureDecisionResults *)
                                                     in_results_wrapper_ptr->object_ptr;
        PictureParentControlSet *pcs = (PictureParentControlSet *)
                                               in_results_ptr->pcs_wrapper->object_ptr;
        SequenceControlSet *scs = pcs->scs;
        if (in_results_ptr->task_type == TASK_TFME)
            me_context_ptr->me_ctx->me_type = ME_MCTF;
        else if (in_results_ptr->task_type == TASK_FIRST_PASS_ME)
            me_context_ptr->me_ctx->me_type = ME_FIRST_PASS;
        else // TASK_PAME or TASK_SUPERRES_RE_ME
            me_context_ptr->me_ctx->me_type = ME_OPEN_LOOP;

        // ME Kernel Signal(s) derivation
        if ((in_results_ptr->task_type == TASK_PAME) ||
            (in_results_ptr->task_type == TASK_SUPERRES_RE_ME))
            if (scs->static_config.pass == ENC_FIRST_PASS)
                first_pass_signal_derivation_me_kernel(scs, pcs, me_context_ptr);
            else
                svt_aom_signal_derivation_me_kernel_oq(scs, pcs, me_context_ptr);

        else if (in_results_ptr->task_type == TASK_TFME)
            tf_signal_derivation_me_kernel_oq(pcs, me_context_ptr);
        else // TASK_FIRST_PASS_ME
            first_pass_signal_derivation_me_kernel(scs, pcs, me_context_ptr);

        if ((in_results_ptr->task_type == TASK_PAME) ||
            (in_results_ptr->task_type == TASK_SUPERRES_RE_ME)) {
            EbPictureBufferDesc *sixteenth_picture_ptr;
            EbPictureBufferDesc *quarter_picture_ptr;
            EbPictureBufferDesc *input_padded_pic;
            EbPictureBufferDesc *input_pic;
            EbPaReferenceObject *pa_ref_obj_;

            //assert((int)pcs->pa_ref_pic_wrapper->live_count > 0);
            pa_ref_obj_ = (EbPaReferenceObject *)
                              pcs->pa_ref_pic_wrapper->object_ptr;
            // Set 1/4 and 1/16 ME input buffer(s); filtered or decimated
            quarter_picture_ptr = (EbPictureBufferDesc *)
                                      pa_ref_obj_->quarter_downsampled_picture_ptr;
            sixteenth_picture_ptr = (EbPictureBufferDesc *)
                                        pa_ref_obj_->sixteenth_downsampled_picture_ptr;
            input_padded_pic = (EbPictureBufferDesc *)pa_ref_obj_->input_padded_pic;

            input_pic = pcs->enhanced_pic;

            // Segments
            uint32_t segment_index   = in_results_ptr->segment_index;
            uint32_t pic_width_in_b64 = (pcs->aligned_width + scs->b64_size - 1) / scs->b64_size;
            uint32_t picture_height_in_b64 = (pcs->aligned_height + scs->b64_size - 1) / scs->b64_size;
            uint32_t y_segment_index;
            uint32_t x_segment_index;

            SEGMENT_CONVERT_IDX_TO_XY(segment_index, x_segment_index, y_segment_index, pcs->me_segments_column_count);
            uint32_t x_b64_start_index = SEGMENT_START_IDX(x_segment_index, pic_width_in_b64, pcs->me_segments_column_count);
            uint32_t x_b64_end_index = SEGMENT_END_IDX(x_segment_index, pic_width_in_b64, pcs->me_segments_column_count);
            uint32_t y_b64_start_index = SEGMENT_START_IDX(y_segment_index, picture_height_in_b64, pcs->me_segments_row_count);
            uint32_t y_b64_end_index = SEGMENT_END_IDX(y_segment_index, picture_height_in_b64, pcs->me_segments_row_count);

            Bool skip_me = FALSE;
            if (scs->static_config.pass == ENC_FIRST_PASS ||
                svt_aom_is_pic_skipped(pcs))
                skip_me = TRUE;
            // skip me for the first pass. ME is already performed
            if (!skip_me) {
                if (pcs->slice_type != I_SLICE) {
                    // Use scaled source references if resolution of the reference is different that of the input
                    svt_aom_use_scaled_source_refs_if_needed(pcs,
                                                     input_pic,
                                                     pa_ref_obj_,
                                                     &input_padded_pic,
                                                     &quarter_picture_ptr,
                                                     &sixteenth_picture_ptr);

                    // 64x64 Block Loop
                    for (uint32_t y_b64_index = y_b64_start_index; y_b64_index < y_b64_end_index; ++y_b64_index) {
                        for (uint32_t x_b64_index = x_b64_start_index; x_b64_index < x_b64_end_index; ++x_b64_index) {

                            uint32_t b64_index    = (uint16_t)(x_b64_index + y_b64_index * pic_width_in_b64);

                            uint32_t b64_origin_x = x_b64_index * scs->b64_size;
                            uint32_t b64_origin_y = y_b64_index * scs->b64_size;

                            // Load the 64x64 Block from the input to the intermediate block buffer
                            uint32_t buffer_index = (input_pic->org_y + b64_origin_y) * input_pic->stride_y +
                                input_pic->org_x + b64_origin_x;
#ifdef ARCH_X86_64
                            uint8_t *src_ptr   = &input_padded_pic->buffer_y[buffer_index];
                            uint32_t b64_height = (pcs->aligned_height - b64_origin_y) < BLOCK_SIZE_64
                                ? pcs->aligned_height - b64_origin_y : BLOCK_SIZE_64;
                            //_MM_HINT_T0     //_MM_HINT_T1    //_MM_HINT_T2//_MM_HINT_NTA
                            for (uint32_t i = 0; i < b64_height; i++) {
                                char const *p = (char const *)(src_ptr + i * input_padded_pic->stride_y);
                                _mm_prefetch(p, _MM_HINT_T2);
                            }
#endif
                            me_context_ptr->me_ctx->b64_src_ptr = &input_padded_pic->buffer_y[buffer_index];
                            me_context_ptr->me_ctx->b64_src_stride = input_padded_pic->stride_y;

                            // Load the 1/4 decimated SB from the 1/4 decimated input to the 1/4 intermediate SB buffer
                            if (me_context_ptr->me_ctx->enable_hme_level1_flag) {
                                buffer_index = (quarter_picture_ptr->org_y + (b64_origin_y >> 1)) * quarter_picture_ptr->stride_y +
                                    quarter_picture_ptr->org_x + (b64_origin_x >> 1);

                                me_context_ptr->me_ctx->quarter_b64_buffer = &quarter_picture_ptr->buffer_y[buffer_index];
                                me_context_ptr->me_ctx->quarter_b64_buffer_stride = quarter_picture_ptr->stride_y;
                            }

                            // Load the 1/16 decimated SB from the 1/16 decimated input to the 1/16 intermediate SB buffer
                            if (me_context_ptr->me_ctx->enable_hme_level0_flag) {
                                buffer_index = (sixteenth_picture_ptr->org_y + (b64_origin_y >> 2)) * sixteenth_picture_ptr->stride_y +
                                    sixteenth_picture_ptr->org_x + (b64_origin_x >> 2);

                                me_context_ptr->me_ctx->sixteenth_b64_buffer = &sixteenth_picture_ptr->buffer_y[buffer_index];
                                me_context_ptr->me_ctx->sixteenth_b64_buffer_stride = sixteenth_picture_ptr->stride_y;
                            }

                            me_context_ptr->me_ctx->me_type = ME_OPEN_LOOP;

                            if ((in_results_ptr->task_type == TASK_PAME) || (in_results_ptr->task_type == TASK_SUPERRES_RE_ME)) {
                                me_context_ptr->me_ctx->num_of_list_to_search =
                                    (pcs->slice_type == P_SLICE) ? 1 /*List 0 only*/
                                    : 2 /*List 0 + 1*/;

                                me_context_ptr->me_ctx->num_of_ref_pic_to_search[0] = pcs->ref_list0_count_try;
                                if (pcs->slice_type == B_SLICE)
                                    me_context_ptr->me_ctx->num_of_ref_pic_to_search[1] = pcs->ref_list1_count_try;
                                me_context_ptr->me_ctx->temporal_layer_index = pcs->temporal_layer_index;
                                me_context_ptr->me_ctx->is_ref = pcs->is_ref;

                                if (pcs->frame_superres_enabled || pcs->frame_resize_enabled) {
                                    for (int i = 0;  i < me_context_ptr->me_ctx->num_of_list_to_search; i++) {
                                        for (int j = 0; j < me_context_ptr->me_ctx->num_of_ref_pic_to_search[i]; j++) {
                                            //assert((int)pcs->ref_pa_pic_ptr_array[i][j]->live_count > 0);
                                            uint8_t sr_denom_idx = svt_aom_get_denom_idx(pcs->superres_denom);
                                            uint8_t resize_denom_idx = svt_aom_get_denom_idx(pcs->resize_denom);
                                            EbPaReferenceObject *ref_object =
                                                (EbPaReferenceObject *)pcs->ref_pa_pic_ptr_array[i][j]->object_ptr;
                                            me_context_ptr->me_ctx->me_ds_ref_array[i][j].picture_ptr =
                                                ref_object->downscaled_input_padded_picture_ptr[sr_denom_idx][resize_denom_idx];
                                            me_context_ptr->me_ctx->me_ds_ref_array[i][j].quarter_picture_ptr =
                                                ref_object->downscaled_quarter_downsampled_picture_ptr[sr_denom_idx][resize_denom_idx];
                                            me_context_ptr->me_ctx->me_ds_ref_array[i][j].sixteenth_picture_ptr =
                                                ref_object->downscaled_sixteenth_downsampled_picture_ptr[sr_denom_idx][resize_denom_idx];
                                            me_context_ptr->me_ctx->me_ds_ref_array[i][j].picture_number =
                                                ref_object->picture_number;
                                        }
                                    }
                                } else {
                                    for (int i = 0; i < me_context_ptr->me_ctx->num_of_list_to_search; i++) {
                                        for (int j = 0; j < me_context_ptr->me_ctx->num_of_ref_pic_to_search[i]; j++) {
                                            //assert((int)pcs->ref_pa_pic_ptr_array[i][j]->live_count > 0);
                                            EbPaReferenceObject *ref_object =
                                                (EbPaReferenceObject *)pcs->ref_pa_pic_ptr_array[i][j]->object_ptr;
                                            me_context_ptr->me_ctx->me_ds_ref_array[i][j].picture_ptr =
                                                ref_object->input_padded_pic;
                                            me_context_ptr->me_ctx->me_ds_ref_array[i][j].quarter_picture_ptr =
                                                ref_object->quarter_downsampled_picture_ptr;
                                            me_context_ptr->me_ctx->me_ds_ref_array[i][j].sixteenth_picture_ptr =
                                                ref_object->sixteenth_downsampled_picture_ptr;
                                            me_context_ptr->me_ctx->me_ds_ref_array[i][j].picture_number =
                                                ref_object->picture_number;
                                        }
                                    }
                                }
                            }

                            svt_aom_motion_estimation_b64(pcs,
                                b64_index,
                                b64_origin_x,
                                b64_origin_y,
                                me_context_ptr->me_ctx,
                                input_pic);

                            if ((in_results_ptr->task_type == TASK_PAME) || (in_results_ptr->task_type == TASK_SUPERRES_RE_ME)) {
                                svt_block_on_mutex(pcs->me_processed_b64_mutex);
                                pcs->me_processed_b64_count++;
                                // We need to finish ME for all SBs to do GM
                                if (pcs->me_processed_b64_count == pcs->b64_total_count) {
                                    if (pcs->gm_ctrls.enabled)
                                        svt_aom_global_motion_estimation(pcs, input_pic);
                                    else
                                        // Initilize global motion to be OFF when GM is OFF
                                        memset(pcs->is_global_motion, FALSE, MAX_NUM_OF_REF_PIC_LIST * REF_LIST_MAX_DEPTH);
                                }

                                svt_release_mutex(pcs->me_processed_b64_mutex);
                            }
                        }
                    }
                }

                if (scs->in_loop_ois == 0 && pcs->tpl_ctrls.enable)
                    for (uint32_t y_b64_index = y_b64_start_index; y_b64_index < y_b64_end_index; ++y_b64_index)
                        for (uint32_t x_b64_index = x_b64_start_index; x_b64_index < x_b64_end_index; ++x_b64_index) {
                            uint32_t b64_index = (uint16_t)(x_b64_index + y_b64_index * pic_width_in_b64);
                            svt_aom_open_loop_intra_search_mb(pcs, b64_index, input_pic);
                        }
            }
            // Get Empty Results Object
            svt_get_empty_object(me_context_ptr->motion_estimation_results_output_fifo_ptr,
                                 &out_results_wrapper);

            MotionEstimationResults *out_results = (MotionEstimationResults *)
                                                           out_results_wrapper->object_ptr;
            out_results->pcs_wrapper = in_results_ptr->pcs_wrapper;
            out_results->segment_index   = segment_index;
            out_results->task_type       = in_results_ptr->task_type;
            // Release the Input Results
            svt_release_object(in_results_wrapper_ptr);

            // Post the Full Results Object
            svt_post_full_object(out_results_wrapper);
        } else if (in_results_ptr->task_type == TASK_TFME) {
            // temporal filtering start
            me_context_ptr->me_ctx->me_type = ME_MCTF;
            svt_av1_init_temporal_filtering(
                pcs->temp_filt_pcs_list, pcs, me_context_ptr, in_results_ptr->segment_index);

            // Release the Input Results
            svt_release_object(in_results_wrapper_ptr);
        } else {
            // first pass start
            me_context_ptr->me_ctx->me_type = ME_FIRST_PASS;
            open_loop_first_pass(pcs, me_context_ptr, in_results_ptr->segment_index);

            // Release the Input Results
            svt_release_object(in_results_wrapper_ptr);
        }
    }

    return NULL;
}
// clang-format on
