#include "EncModeConfig.h"
#include <stdlib.h>

#include "EbRateDistortionCost.h"
#include "aom_dsp_rtcd.h"
#include "EbCodingLoop.h"

#define LOW_8x8_DIST_VAR_TH 25000
#define HIGH_8x8_DIST_VAR_TH 50000

static uint8_t pf_gi[16] = {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60};

// use this function to set the enable_me_8x8 level
uint8_t svt_aom_get_enable_me_8x8(EncMode enc_mode, bool rtc_tune) {
    uint8_t enable_me_8x8 = 0;
    if ((!rtc_tune && enc_mode <= ENC_M11) || (rtc_tune && enc_mode <= ENC_M10))
        enable_me_8x8 = 1;
    else
        enable_me_8x8 = 0;

    return enable_me_8x8;
}
uint8_t svt_aom_get_enable_me_16x16(EncMode enc_mode, bool rtc_tune) {
    uint8_t enable_me_16x16;
    if ((enc_mode <= ENC_M13) || (rtc_tune && (enc_mode <= ENC_M13)))
        enable_me_16x16 = 1;
    else
        enable_me_16x16 = 0;

    return enable_me_16x16;
}

uint8_t svt_aom_get_gm_core_level(EncMode enc_mode, bool super_res_off) {
    uint8_t gm_level = 0;
    if (super_res_off) {
        if (enc_mode <= ENC_M1)
            gm_level = 2;
        else if (enc_mode <= ENC_M2)
            gm_level = 3;
        else if (enc_mode <= ENC_M4)
            gm_level = 4;
        else if (enc_mode <= ENC_M6)
            gm_level = 5;
        else
            gm_level = 0;
    }
    return gm_level;
}
bool svt_aom_need_gm_ref_info(EncMode enc_mode, bool super_res_off) {
    uint8_t                 gm_lvl = svt_aom_get_gm_core_level(enc_mode, super_res_off);
    PictureParentControlSet pcs_tmp;
    svt_aom_set_gm_controls(&pcs_tmp, gm_lvl);
    return pcs_tmp.gm_ctrls.use_ref_info;
}

uint8_t svt_aom_derive_gm_level(PictureParentControlSet *pcs, bool super_res_off) {
    SequenceControlSet *scs       = pcs->scs;
    uint8_t             gm_level  = 0;
    const EncMode       enc_mode  = pcs->enc_mode;
    const uint8_t       is_islice = pcs->slice_type == I_SLICE;
    // disable global motion when reference scaling enabled,
    // even if current pic is not scaled, because its reference
    // pics might be scaled in different size
    // super-res is ok for its reference pics are always upscaled
    // to original size
    if (scs->enable_global_motion && !is_islice)
        gm_level = svt_aom_get_gm_core_level(enc_mode, super_res_off);
    return gm_level;
}
/************************************************
 * Set ME/HME Params from Config
 ************************************************/
/************************************************
 * Set HME Search area parameters
 ************************************************/
static void set_hme_search_params(PictureParentControlSet *pcs, MeContext *me_ctx, EbInputResolution input_resolution) {
    const bool rtc_tune = (pcs->scs->static_config.pred_structure == SVT_AV1_PRED_LOW_DELAY_B) ? true : false;
    // Set number of HME level 0 search regions to use
    me_ctx->num_hme_sa_w = 2;
    me_ctx->num_hme_sa_h = 2;

    // Set HME level 0 min and max search areas
    if (pcs->enc_mode <= ENC_MR) {
        if (input_resolution < INPUT_SIZE_4K_RANGE) {
            me_ctx->hme_l0_sa.sa_min = (SearchArea){32, 32};
            me_ctx->hme_l0_sa.sa_max = (SearchArea){192, 192};
        } else {
            me_ctx->hme_l0_sa.sa_min = (SearchArea){240, 240};
            me_ctx->hme_l0_sa.sa_max = (SearchArea){480, 480};
        }
    } else if (pcs->enc_mode <= ENC_M4) {
        me_ctx->hme_l0_sa.sa_min = (SearchArea){32, 32};
        me_ctx->hme_l0_sa.sa_max = (SearchArea){192, 192};
    } else if (pcs->enc_mode <= ENC_M6) {
        if (input_resolution <= INPUT_SIZE_1080p_RANGE) {
            if (pcs->sc_class1) {
                me_ctx->hme_l0_sa.sa_min = (SearchArea){32, 32};
                me_ctx->hme_l0_sa.sa_max = (SearchArea){192, 192};
            } else {
                me_ctx->hme_l0_sa.sa_min = (SearchArea){16, 16};
                me_ctx->hme_l0_sa.sa_max = (SearchArea){192, 192};
            }
        } else {
            me_ctx->hme_l0_sa.sa_min = (SearchArea){32, 32};
            me_ctx->hme_l0_sa.sa_max = (SearchArea){192, 192};
        }
    } else if ((!rtc_tune && pcs->enc_mode <= ENC_M10) || (rtc_tune && pcs->enc_mode <= ENC_M9)) {
        if (pcs->sc_class1) {
            me_ctx->hme_l0_sa.sa_min = (SearchArea){32, 32};
            me_ctx->hme_l0_sa.sa_max = (SearchArea){192, 192};
        } else {
            me_ctx->hme_l0_sa.sa_min = (SearchArea){16, 16};
            me_ctx->hme_l0_sa.sa_max = (SearchArea){192, 192};
        }
    } else {
        if (pcs->sc_class1) {
            me_ctx->hme_l0_sa.sa_min = (SearchArea){32, 32};
            me_ctx->hme_l0_sa.sa_max = (SearchArea){192, 192};
        } else {
            if (input_resolution < INPUT_SIZE_4K_RANGE) {
                me_ctx->hme_l0_sa.sa_min = (SearchArea){8, 8};
                me_ctx->hme_l0_sa.sa_max = (SearchArea){96, 96};
            } else {
                me_ctx->hme_l0_sa.sa_min = (SearchArea){16, 16};
                me_ctx->hme_l0_sa.sa_max = (SearchArea){96, 96};
            }
        }
    }
    // Set the HME Level 1 and Level 2 refinement areas
    if (pcs->enc_mode <= ENC_M0) {
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
static void set_me_search_params(SequenceControlSet *scs, PictureParentControlSet *pcs, MeContext *me_ctx,
                                 EbInputResolution input_resolution) {
    uint32_t   hierarchical_levels = pcs->hierarchical_levels;
    const bool rtc_tune            = (scs->static_config.pred_structure == SVT_AV1_PRED_LOW_DELAY_B) ? true : false;
    // Set the min and max ME search area
    if (rtc_tune) {
        if (pcs->sc_class1) {
            if (pcs->enc_mode <= ENC_M9) {
                me_ctx->me_sa.sa_min = (SearchArea){32, 32};
                me_ctx->me_sa.sa_max = (SearchArea){96, 96};
            } else if (pcs->enc_mode <= ENC_M11) {
                if (input_resolution < INPUT_SIZE_1080p_RANGE) {
                    me_ctx->me_sa.sa_min = (SearchArea){16, 16};
                    me_ctx->me_sa.sa_max = (SearchArea){32, 16};
                } else {
                    me_ctx->me_sa.sa_min = (SearchArea){24, 24};
                    me_ctx->me_sa.sa_max = (SearchArea){24, 24};
                }
            } else {
                if (input_resolution < INPUT_SIZE_1080p_RANGE) {
                    me_ctx->me_sa.sa_min = (SearchArea){16, 16};
                    me_ctx->me_sa.sa_max = (SearchArea){32, 16};
                } else {
                    me_ctx->me_sa.sa_min = (SearchArea){16, 6};
                    me_ctx->me_sa.sa_max = (SearchArea){16, 9};
                }
            }
        } else {
            if (pcs->enc_mode <= ENC_M10) {
                if (input_resolution < INPUT_SIZE_1080p_RANGE) {
                    me_ctx->me_sa.sa_min = (SearchArea){16, 16};
                    me_ctx->me_sa.sa_max = (SearchArea){32, 16};
                } else {
                    me_ctx->me_sa.sa_min = (SearchArea){16, 6};
                    me_ctx->me_sa.sa_max = (SearchArea){16, 9};
                }
            } else {
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
            }
        }
    } else if (pcs->sc_class1) {
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
            me_ctx->me_sa.sa_min = (SearchArea){16, 16};
            me_ctx->me_sa.sa_max = (SearchArea){96, 96};
        } else {
            me_ctx->me_sa.sa_min = (SearchArea){8, 8};
            me_ctx->me_sa.sa_max = (SearchArea){32, 32};
        }
    } else if (pcs->enc_mode <= ENC_MR) {
        me_ctx->me_sa.sa_min = (SearchArea){64, 64};
        me_ctx->me_sa.sa_max = (SearchArea){256, 256};
    } else if (pcs->enc_mode <= ENC_M2) {
        me_ctx->me_sa.sa_min = (SearchArea){64, 64};
        me_ctx->me_sa.sa_max = (SearchArea){128, 128};
    } else if (pcs->enc_mode <= ENC_M3) {
        me_ctx->me_sa.sa_min = (SearchArea){16, 16};
        me_ctx->me_sa.sa_max = (SearchArea){128, 128};
    } else if (pcs->enc_mode <= ENC_M7) {
        me_ctx->me_sa.sa_min = (SearchArea){16, 16};
        me_ctx->me_sa.sa_max = (SearchArea){64, 32};
    } else if (pcs->enc_mode <= ENC_M8) {
        if (hierarchical_levels <= 3) {
            if (input_resolution < INPUT_SIZE_4K_RANGE) {
                me_ctx->me_sa.sa_min = (SearchArea){8, 5};
                me_ctx->me_sa.sa_max = (SearchArea){16, 9};
            } else {
                me_ctx->me_sa.sa_min = (SearchArea){8, 1};
                me_ctx->me_sa.sa_max = (SearchArea){8, 1};
            }
        } else if (input_resolution < INPUT_SIZE_1080p_RANGE) {
            me_ctx->me_sa.sa_min = (SearchArea){16, 16};
            me_ctx->me_sa.sa_max = (SearchArea){32, 16};
        } else {
            me_ctx->me_sa.sa_min = (SearchArea){16, 6};
            me_ctx->me_sa.sa_max = (SearchArea){16, 9};
        }
    } else if (pcs->enc_mode <= ENC_M10) {
        if (input_resolution < INPUT_SIZE_1080p_RANGE) {
            me_ctx->me_sa.sa_min = (SearchArea){16, 16};
            me_ctx->me_sa.sa_max = (SearchArea){32, 16};
        } else {
            me_ctx->me_sa.sa_min = (SearchArea){16, 6};
            me_ctx->me_sa.sa_max = (SearchArea){16, 9};
        }
    } else if (pcs->enc_mode <= ENC_M12) {
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
    // Scale up the MIN ME area if low frame rate
    bool low_frame_rate_flag = (scs->frame_rate >> 16);
    if (low_frame_rate_flag) {
        me_ctx->me_sa.sa_min.width  = (me_ctx->me_sa.sa_min.width * 3) >> 1;
        me_ctx->me_sa.sa_min.height = (me_ctx->me_sa.sa_min.height * 3) >> 1;
    }
}
static void svt_aom_set_me_hme_ref_prune_ctrls(MeContext *me_ctx, uint8_t prune_level) {
    MeHmeRefPruneCtrls *me_hme_prune_ctrls = &me_ctx->me_hme_prune_ctrls;

    switch (prune_level) {
    case 0:
        me_hme_prune_ctrls->enable_me_hme_ref_pruning               = 0;
        me_hme_prune_ctrls->prune_ref_if_hme_sad_dev_bigger_than_th = (uint16_t)~0;
        me_hme_prune_ctrls->prune_ref_if_me_sad_dev_bigger_than_th  = (uint16_t)~0;

        me_hme_prune_ctrls->zz_sad_th    = 0;
        me_hme_prune_ctrls->zz_sad_pct   = 0;
        me_hme_prune_ctrls->phme_sad_th  = 0;
        me_hme_prune_ctrls->phme_sad_pct = 0;
        break;
    case 1:
        me_hme_prune_ctrls->enable_me_hme_ref_pruning               = 1;
        me_hme_prune_ctrls->prune_ref_if_hme_sad_dev_bigger_than_th = 80;
        me_hme_prune_ctrls->prune_ref_if_me_sad_dev_bigger_than_th  = 60;
        me_hme_prune_ctrls->protect_closest_refs                    = 1;

        me_hme_prune_ctrls->zz_sad_th    = 0;
        me_hme_prune_ctrls->zz_sad_pct   = 0;
        me_hme_prune_ctrls->phme_sad_th  = 0;
        me_hme_prune_ctrls->phme_sad_pct = 0;
        break;
    case 2:
        me_hme_prune_ctrls->enable_me_hme_ref_pruning               = 1;
        me_hme_prune_ctrls->prune_ref_if_hme_sad_dev_bigger_than_th = 5;
        me_hme_prune_ctrls->prune_ref_if_me_sad_dev_bigger_than_th  = 60;
        me_hme_prune_ctrls->protect_closest_refs                    = 1;

        me_hme_prune_ctrls->zz_sad_th    = 0;
        me_hme_prune_ctrls->zz_sad_pct   = 0;
        me_hme_prune_ctrls->phme_sad_th  = 0;
        me_hme_prune_ctrls->phme_sad_pct = 0;
        break;
    case 3:
        me_hme_prune_ctrls->enable_me_hme_ref_pruning               = 1;
        me_hme_prune_ctrls->prune_ref_if_hme_sad_dev_bigger_than_th = 5;
        me_hme_prune_ctrls->prune_ref_if_me_sad_dev_bigger_than_th  = 60;
        me_hme_prune_ctrls->protect_closest_refs                    = 1;

        me_hme_prune_ctrls->zz_sad_th    = 20 * 64 * 64;
        me_hme_prune_ctrls->zz_sad_pct   = 5;
        me_hme_prune_ctrls->phme_sad_th  = 0;
        me_hme_prune_ctrls->phme_sad_pct = 0;
        break;
    case 4:
        me_hme_prune_ctrls->enable_me_hme_ref_pruning               = 1;
        me_hme_prune_ctrls->prune_ref_if_hme_sad_dev_bigger_than_th = 5;
        me_hme_prune_ctrls->prune_ref_if_me_sad_dev_bigger_than_th  = 60;
        me_hme_prune_ctrls->protect_closest_refs                    = 1;

        me_hme_prune_ctrls->zz_sad_th    = 20 * 64 * 64;
        me_hme_prune_ctrls->zz_sad_pct   = 5;
        me_hme_prune_ctrls->phme_sad_th  = 10 * 64 * 64;
        me_hme_prune_ctrls->phme_sad_pct = 5;
        break;
    case 5:
        me_hme_prune_ctrls->enable_me_hme_ref_pruning               = 1;
        me_hme_prune_ctrls->prune_ref_if_hme_sad_dev_bigger_than_th = 200;
        me_hme_prune_ctrls->prune_ref_if_me_sad_dev_bigger_than_th  = 60;
        me_hme_prune_ctrls->protect_closest_refs                    = 0;

        me_hme_prune_ctrls->zz_sad_th    = 20 * 64 * 64;
        me_hme_prune_ctrls->zz_sad_pct   = 5;
        me_hme_prune_ctrls->phme_sad_th  = 10 * 64 * 64;
        me_hme_prune_ctrls->phme_sad_pct = 5;
        break;
    default: assert(0); break;
    }
}

static void svt_aom_set_me_sr_adjustment_ctrls(MeContext *me_ctx, uint8_t sr_adjustment_level) {
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
        me_sr_adjustment_ctrls->enable_me_sr_adjustment              = 2;
        me_sr_adjustment_ctrls->reduce_me_sr_based_on_mv_length_th   = 16;
        me_sr_adjustment_ctrls->stationary_hme_sad_abs_th            = 20000;
        me_sr_adjustment_ctrls->stationary_me_sr_divisor             = 8;
        me_sr_adjustment_ctrls->reduce_me_sr_based_on_hme_sad_abs_th = 20000;
        me_sr_adjustment_ctrls->me_sr_divisor_for_low_hme_sad        = 8;
        me_sr_adjustment_ctrls->distance_based_hme_resizing          = 1;
        break;

    case 5:

        me_sr_adjustment_ctrls->enable_me_sr_adjustment              = 2;
        me_sr_adjustment_ctrls->reduce_me_sr_based_on_mv_length_th   = 20;
        me_sr_adjustment_ctrls->stationary_hme_sad_abs_th            = 24000;
        me_sr_adjustment_ctrls->stationary_me_sr_divisor             = 8;
        me_sr_adjustment_ctrls->reduce_me_sr_based_on_hme_sad_abs_th = 24000;
        me_sr_adjustment_ctrls->me_sr_divisor_for_low_hme_sad        = 8;
        me_sr_adjustment_ctrls->distance_based_hme_resizing          = 1;

        break;

    default: assert(0); break;
    }
    if (me_ctx->enable_hme_level2_flag == 0) {
        if (me_ctx->enable_hme_level1_flag == 1) {
            me_sr_adjustment_ctrls->stationary_hme_sad_abs_th = me_sr_adjustment_ctrls->stationary_hme_sad_abs_th / 4;
            me_sr_adjustment_ctrls->reduce_me_sr_based_on_hme_sad_abs_th =
                me_sr_adjustment_ctrls->reduce_me_sr_based_on_hme_sad_abs_th / 4;
        } else {
            me_sr_adjustment_ctrls->stationary_hme_sad_abs_th = me_sr_adjustment_ctrls->stationary_hme_sad_abs_th / 16;
            me_sr_adjustment_ctrls->reduce_me_sr_based_on_hme_sad_abs_th =
                me_sr_adjustment_ctrls->reduce_me_sr_based_on_hme_sad_abs_th / 16;
        }
    }
}
static void svt_aom_set_me_8x8_var_ctrls(MeContext *me_ctx, uint8_t level) {
    Me8x8VarCtrls *me_8x8_var_ctrls = &me_ctx->me_8x8_var_ctrls;

    switch (level) {
    case 0: me_8x8_var_ctrls->enabled = 0; break;
    case 1:
        me_8x8_var_ctrls->enabled       = 1;
        me_8x8_var_ctrls->me_sr_div4_th = 80000;
        me_8x8_var_ctrls->me_sr_div2_th = 150000;
        break;
    default: assert(0);
    }
}
/*configure PreHme control*/
static void svt_aom_set_prehme_ctrls(MeContext *me_ctx, uint8_t level) {
    PreHmeCtrls *ctrl = &me_ctx->prehme_ctrl;

    switch (level) {
    case 0: ctrl->enable = 0; break;
    case 1:
        ctrl->enable = 1;
        // vertical shape search region
        ctrl->prehme_sa_cfg[0].sa_min = (SearchArea){8, 144};
        ctrl->prehme_sa_cfg[0].sa_max = (SearchArea){8, 496};
        // horizontal shape search region
        ctrl->prehme_sa_cfg[1].sa_min = (SearchArea){144, 3};
        ctrl->prehme_sa_cfg[1].sa_max = (SearchArea){496, 3};
        ctrl->skip_search_line        = 0;
        ctrl->l1_early_exit           = 0;
        break;
    case 2:
        ctrl->enable = 1;
        // vertical shape search region
        ctrl->prehme_sa_cfg[0].sa_min = (SearchArea){8, 100};
        ctrl->prehme_sa_cfg[0].sa_max = (SearchArea){8, 400};
        // horizontal shape search region
        ctrl->prehme_sa_cfg[1].sa_min = (SearchArea){96, 3};
        ctrl->prehme_sa_cfg[1].sa_max = (SearchArea){384, 3};
        ctrl->skip_search_line        = 0;
        ctrl->l1_early_exit           = 0;
        break;
    case 3:
        ctrl->enable = 1;
        // vertical shape search region
        ctrl->prehme_sa_cfg[0].sa_min = (SearchArea){8, 100};
        ctrl->prehme_sa_cfg[0].sa_max = (SearchArea){8, 350};
        // horizontal shape search region
        ctrl->prehme_sa_cfg[1].sa_min = (SearchArea){32, 7};
        ctrl->prehme_sa_cfg[1].sa_max = (SearchArea){200, 7};
        ctrl->skip_search_line        = 1;
        ctrl->l1_early_exit           = 0;
        break;
    case 4:
        ctrl->enable = 1;
        // vertical shape search region
        ctrl->prehme_sa_cfg[0].sa_min = (SearchArea){8, 100};
        ctrl->prehme_sa_cfg[0].sa_max = (SearchArea){8, 350};
        // horizontal shape search region
        ctrl->prehme_sa_cfg[1].sa_min = (SearchArea){32, 7};
        ctrl->prehme_sa_cfg[1].sa_max = (SearchArea){128, 7};
        ctrl->skip_search_line        = 1;
        ctrl->l1_early_exit           = 1;
        break;
    default: assert(0); break;
    }
}

/************************************************
 * Set ME/HME Params for Altref Temporal Filtering
 ************************************************/
static void tf_set_me_hme_params_oq(MeContext *me_ctx, PictureParentControlSet *pcs) {
    switch (pcs->tf_ctrls.hme_me_level) {
    case 0:
        me_ctx->num_hme_sa_w                = 2;
        me_ctx->num_hme_sa_h                = 2;
        me_ctx->hme_l0_sa_default_tf.sa_min = (SearchArea){30, 30};
        me_ctx->hme_l0_sa_default_tf.sa_max = (SearchArea){60, 60};
        me_ctx->hme_l1_sa                   = (SearchArea){16, 16};
        me_ctx->hme_l2_sa                   = (SearchArea){16, 16};
        me_ctx->me_sa.sa_min                = (SearchArea){60, 60};
        me_ctx->me_sa.sa_max                = (SearchArea){120, 120};
        break;

    case 1:
        me_ctx->num_hme_sa_w                = 2;
        me_ctx->num_hme_sa_h                = 2;
        me_ctx->hme_l0_sa_default_tf.sa_min = (SearchArea){16, 16};
        me_ctx->hme_l0_sa_default_tf.sa_max = (SearchArea){32, 32};
        me_ctx->hme_l1_sa                   = (SearchArea){16, 16};
        me_ctx->hme_l2_sa                   = (SearchArea){16, 16};
        me_ctx->me_sa.sa_min                = (SearchArea){16, 16};
        me_ctx->me_sa.sa_max                = (SearchArea){32, 32};
        break;

    case 2:
        me_ctx->num_hme_sa_w                = 2;
        me_ctx->num_hme_sa_h                = 2;
        me_ctx->hme_l0_sa_default_tf.sa_min = (SearchArea){8, 8};
        me_ctx->hme_l0_sa_default_tf.sa_max = (SearchArea){16, 16};
        me_ctx->hme_l1_sa                   = (SearchArea){16, 16};
        me_ctx->hme_l2_sa                   = (SearchArea){16, 16};
        me_ctx->me_sa.sa_min                = (SearchArea){8, 4};
        me_ctx->me_sa.sa_max                = (SearchArea){16, 8};
        break;
    case 3:
        me_ctx->num_hme_sa_w                = 2;
        me_ctx->num_hme_sa_h                = 2;
        me_ctx->hme_l0_sa_default_tf.sa_min = (SearchArea){4, 4};
        me_ctx->hme_l0_sa_default_tf.sa_max = (SearchArea){4, 4};
        me_ctx->hme_l1_sa                   = (SearchArea){8, 8};
        me_ctx->hme_l2_sa                   = (SearchArea){8, 8};
        me_ctx->me_sa.sa_min                = (SearchArea){8, 8};
        me_ctx->me_sa.sa_max                = (SearchArea){8, 8};
        break;

    default: assert(0); break;
    }
};
/******************************************************
* Derive ME Settings for OQ
  Input   : encoder mode and tune
  Output  : ME Kernel signal(s)
******************************************************/
void svt_aom_sig_deriv_me(SequenceControlSet *scs, PictureParentControlSet *pcs, MeContext *me_ctx) {
    EncMode           enc_mode         = pcs->enc_mode;
    EbInputResolution input_resolution = scs->input_resolution;
    const bool        rtc_tune         = (scs->static_config.pred_structure == SVT_AV1_PRED_LOW_DELAY_B) ? true : false;
    // Set ME search area
    set_me_search_params(scs, pcs, me_ctx, input_resolution);

    // Set HME search area
    set_hme_search_params(pcs, me_ctx, input_resolution);

    // Set HME flags
    me_ctx->enable_hme_flag        = pcs->enable_hme_flag;
    me_ctx->enable_hme_level0_flag = pcs->enable_hme_level0_flag;
    me_ctx->enable_hme_level1_flag = pcs->enable_hme_level1_flag;
    me_ctx->enable_hme_level2_flag = pcs->enable_hme_level2_flag;
    // HME Search Method
    me_ctx->hme_search_method = SUB_SAD_SEARCH;
    me_ctx->me_search_method  = SUB_SAD_SEARCH;

    if (pcs->sc_class1) {
        if (rtc_tune) {
            me_ctx->reduce_hme_l0_sr_th_min = 8;
            me_ctx->reduce_hme_l0_sr_th_max = 100;
        } else {
            me_ctx->reduce_hme_l0_sr_th_min = 0;
            me_ctx->reduce_hme_l0_sr_th_max = 0;
        }
    } else if ((rtc_tune && pcs->enc_mode <= ENC_M9) || (!rtc_tune && pcs->enc_mode <= ENC_M11)) {
        me_ctx->reduce_hme_l0_sr_th_min = 0;
        me_ctx->reduce_hme_l0_sr_th_max = 0;
    } else {
        me_ctx->reduce_hme_l0_sr_th_min = 8;
        me_ctx->reduce_hme_l0_sr_th_max = 200;
    }
    // Set pre-hme level (0-2)
    uint8_t prehme_level = 0;
    if (pcs->sc_class1)
        prehme_level = rtc_tune ? 1 : 2;
    else if (rtc_tune) {
        if (enc_mode <= ENC_M9)
            prehme_level = 2;
        else if (enc_mode <= ENC_M12)
            prehme_level = 4;
        else
            prehme_level = 0;
    } else {
        if (enc_mode <= ENC_M9)
            prehme_level = 2;
        else
            prehme_level = 4;
    }
    if (pcs->enable_hme_level1_flag == 0)
        prehme_level = 0;

    svt_aom_set_prehme_ctrls(me_ctx, prehme_level);

    // Set hme/me based reference pruning level (0-4)
    if (pcs->sc_class1) {
        if (enc_mode <= ENC_M1)
            svt_aom_set_me_hme_ref_prune_ctrls(me_ctx, 1);
        else if (enc_mode <= ENC_M9)
            svt_aom_set_me_hme_ref_prune_ctrls(me_ctx, 2);
        else
            svt_aom_set_me_hme_ref_prune_ctrls(me_ctx, 5);
    } else {
        if (enc_mode <= ENC_M0)
            svt_aom_set_me_hme_ref_prune_ctrls(me_ctx, 1);
        else if (enc_mode <= ENC_M7) {
            if (pcs->temporal_layer_index == 0)
                svt_aom_set_me_hme_ref_prune_ctrls(me_ctx, 1);
            else
                svt_aom_set_me_hme_ref_prune_ctrls(me_ctx, 2);
        } else if (enc_mode <= ENC_M9) {
            if (pcs->temporal_layer_index == 0)
                svt_aom_set_me_hme_ref_prune_ctrls(me_ctx, 1);
            else
                svt_aom_set_me_hme_ref_prune_ctrls(me_ctx, 3);
        } else if (enc_mode <= ENC_M11)
            svt_aom_set_me_hme_ref_prune_ctrls(me_ctx, 3);
        else
            svt_aom_set_me_hme_ref_prune_ctrls(me_ctx, 4);
    }
    // Set hme-based me sr adjustment level
    if (pcs->sc_class1)
        if ((!rtc_tune && enc_mode <= ENC_M9) || (rtc_tune && enc_mode <= ENC_M8))
            svt_aom_set_me_sr_adjustment_ctrls(me_ctx, 4);
        else
            svt_aom_set_me_sr_adjustment_ctrls(me_ctx, 5);
    else if (enc_mode <= ENC_M2)
        svt_aom_set_me_sr_adjustment_ctrls(me_ctx, 1);
    else
        svt_aom_set_me_sr_adjustment_ctrls(me_ctx, 3);
    uint8_t me_8x8_var_lvl = 0;
    if (enc_mode <= ENC_M0)
        me_8x8_var_lvl = 0;
    else
        me_8x8_var_lvl = 1;
    svt_aom_set_me_8x8_var_ctrls(me_ctx, me_8x8_var_lvl);
    if (enc_mode <= ENC_M2)
        me_ctx->prune_me_candidates_th = 0;
    else
        me_ctx->prune_me_candidates_th = 65;
    // Set signal at picture level b/c may check signal in MD
    me_ctx->use_best_unipred_cand_only = pcs->use_best_me_unipred_cand_only;
    if (rtc_tune) {
        if (pcs->sc_class1)
            me_ctx->me_early_exit_th = BLOCK_SIZE_64 * BLOCK_SIZE_64;
        else if (enc_mode <= ENC_M11)
            me_ctx->me_early_exit_th = BLOCK_SIZE_64 * BLOCK_SIZE_64 * 8;
        else
            me_ctx->me_early_exit_th = BLOCK_SIZE_64 * BLOCK_SIZE_64 * 9;
    } else {
        if (enc_mode <= ENC_M3)
            me_ctx->me_early_exit_th = 0;
        else
            me_ctx->me_early_exit_th = BLOCK_SIZE_64 * BLOCK_SIZE_64 * 8;
    }

    me_ctx->me_safe_limit_zz_th = scs->mrp_ctrls.safe_limit_nref == 1 ? scs->mrp_ctrls.safe_limit_zz_th : 0;

    me_ctx->skip_frame                  = 0;
    me_ctx->prev_me_stage_based_exit_th = 0;
    if (rtc_tune && pcs->sc_class1) {
        me_ctx->prev_me_stage_based_exit_th = BLOCK_SIZE_64 * BLOCK_SIZE_64 * 4;
    } else {
        me_ctx->prev_me_stage_based_exit_th = 0;
    }
};
/******************************************************
* Derive ME Settings for OQ for Altref Temporal Filtering
  Input   : encoder mode and tune
  Output  : ME Kernel signal(s)
******************************************************/
void svt_aom_sig_deriv_me_tf(PictureParentControlSet *pcs, MeContext *me_ctx) {
    const EbInputResolution resolution = pcs->scs->input_resolution;
    const uint8_t           enc_mode   = pcs->enc_mode;
    // Set ME/HME search regions
    tf_set_me_hme_params_oq(me_ctx, pcs);
    // Set HME flags
    me_ctx->enable_hme_flag        = pcs->tf_enable_hme_flag;
    me_ctx->enable_hme_level0_flag = pcs->tf_enable_hme_level0_flag;
    me_ctx->enable_hme_level1_flag = pcs->tf_enable_hme_level1_flag;
    me_ctx->enable_hme_level2_flag = pcs->tf_enable_hme_level2_flag;
    if (pcs->tf_ctrls.hme_me_level <= 1) {
        // HME Search Method
        me_ctx->hme_search_method = FULL_SAD_SEARCH;
        // ME Search Method
        me_ctx->me_search_method = FULL_SAD_SEARCH;
    } else {
        // HME Search Method
        me_ctx->hme_search_method = SUB_SAD_SEARCH;
        // ME Search Method
        me_ctx->me_search_method = SUB_SAD_SEARCH;
    }
    uint8_t prehme_level = 0;
    svt_aom_set_prehme_ctrls(me_ctx, prehme_level);

    // Set hme/me based reference pruning level (0-4)
    // Ref pruning is disallowed for TF in motion_estimate_sb()
    svt_aom_set_me_hme_ref_prune_ctrls(me_ctx, 0);

    // Set hme-based me sr adjustment level
    // ME SR adjustment is disallowed for TF in motion_estimate_sb()
    svt_aom_set_me_sr_adjustment_ctrls(me_ctx, 0);

    svt_aom_set_me_8x8_var_ctrls(me_ctx, 0);
    me_ctx->me_early_exit_th            = enc_mode <= ENC_M6 || resolution <= INPUT_SIZE_720p_RANGE
                   ? 0
                   : BLOCK_SIZE_64 * BLOCK_SIZE_64 * 4;
    me_ctx->me_safe_limit_zz_th         = 0;
    me_ctx->reduce_hme_l0_sr_th_min     = 0;
    me_ctx->reduce_hme_l0_sr_th_max     = 0;
    me_ctx->skip_frame                  = 0;
    me_ctx->prev_me_stage_based_exit_th = enc_mode <= ENC_M6 || resolution <= INPUT_SIZE_720p_RANGE
        ? 0
        : BLOCK_SIZE_64 * BLOCK_SIZE_64 * 4;
};

/* Wrapper function to compute TPL Synthesizer block size: Used in init memory allocation and TPL Controls*/
uint8_t svt_aom_get_tpl_synthesizer_block_size(int8_t tpl_level, uint32_t picture_width, uint32_t picture_height) {
    uint8_t blk_size;
    if (tpl_level <= 6)
        blk_size = 16;
    else
        blk_size = AOMMIN(picture_width, picture_height) >= 720 ? 32 : 16;

    return blk_size;
}

/*************************************************************************************
Set the TPL controls that control TPL search complexity.
***************************************************************************************/
void svt_aom_set_tpl_extended_controls(PictureParentControlSet *pcs, uint8_t tpl_level) {
    TplControls            *tpl_ctrls  = &pcs->tpl_ctrls;
    SequenceControlSet     *scs        = pcs->scs;
    const uint8_t           is_islice  = pcs->slice_type == I_SLICE;
    const EbInputResolution resolution = scs->input_resolution;

    switch (tpl_level) {
    case 0:
        tpl_ctrls->enable                  = 0;
        tpl_ctrls->compute_rate            = 0;
        tpl_ctrls->enable_tpl_qps          = 0;
        tpl_ctrls->disable_intra_pred_nref = 0;
        tpl_ctrls->intra_mode_end          = DC_PRED;
        tpl_ctrls->reduced_tpl_group       = -1;
        tpl_ctrls->pf_shape                = DEFAULT_SHAPE;
        tpl_ctrls->use_sad_in_src_search   = 0;
        tpl_ctrls->dispenser_search_level  = 0;
        tpl_ctrls->subsample_tx            = 0;
        tpl_ctrls->subpel_depth            = FULL_PEL;
        break;
    case 1:
        tpl_ctrls->enable                  = 1;
        tpl_ctrls->compute_rate            = 1;
        tpl_ctrls->enable_tpl_qps          = 1;
        tpl_ctrls->disable_intra_pred_nref = 0;
        tpl_ctrls->intra_mode_end          = PAETH_PRED;
        tpl_ctrls->reduced_tpl_group       = -1;
        tpl_ctrls->pf_shape                = DEFAULT_SHAPE;
        tpl_ctrls->use_sad_in_src_search   = 0;
        tpl_ctrls->dispenser_search_level  = 0;
        tpl_ctrls->subsample_tx            = 0;
        tpl_ctrls->subpel_depth            = QUARTER_PEL;
        break;
    case 2:
        tpl_ctrls->enable                  = 1;
        tpl_ctrls->compute_rate            = 0;
        tpl_ctrls->enable_tpl_qps          = 0;
        tpl_ctrls->disable_intra_pred_nref = 0;
        tpl_ctrls->intra_mode_end          = PAETH_PRED;
        tpl_ctrls->reduced_tpl_group       = is_islice ? -1 : (pcs->hierarchical_levels == 5 ? 4 : 3);
        tpl_ctrls->pf_shape                = resolution <= INPUT_SIZE_480p_RANGE ? N2_SHAPE : N4_SHAPE;
        tpl_ctrls->use_sad_in_src_search   = 1;
        tpl_ctrls->dispenser_search_level  = 0;
        tpl_ctrls->subsample_tx            = 0;
        tpl_ctrls->subpel_depth            = QUARTER_PEL;
        break;
    case 3:
        tpl_ctrls->enable                  = 1;
        tpl_ctrls->compute_rate            = 0;
        tpl_ctrls->enable_tpl_qps          = 0;
        tpl_ctrls->disable_intra_pred_nref = 1;
        tpl_ctrls->intra_mode_end          = DC_PRED;
        tpl_ctrls->reduced_tpl_group       = is_islice ? -1 : (pcs->hierarchical_levels == 5 ? 4 : 3);
        tpl_ctrls->pf_shape                = resolution <= INPUT_SIZE_480p_RANGE ? N2_SHAPE : N4_SHAPE;
        tpl_ctrls->use_sad_in_src_search   = 1;
        tpl_ctrls->dispenser_search_level  = 0;
        tpl_ctrls->subsample_tx            = 0;
        tpl_ctrls->subpel_depth            = FULL_PEL;
        break;
    case 4:
        tpl_ctrls->enable                  = 1;
        tpl_ctrls->compute_rate            = 0;
        tpl_ctrls->enable_tpl_qps          = 0;
        tpl_ctrls->disable_intra_pred_nref = 1;
        tpl_ctrls->intra_mode_end          = DC_PRED;
        tpl_ctrls->reduced_tpl_group       = pcs->hierarchical_levels == 5 ? (is_islice ? 4 : 3) : (is_islice ? 3 : 2);
        tpl_ctrls->pf_shape                = resolution <= INPUT_SIZE_480p_RANGE ? N2_SHAPE : N4_SHAPE;
        tpl_ctrls->use_sad_in_src_search   = 1;
        tpl_ctrls->dispenser_search_level  = 0;
        tpl_ctrls->subsample_tx            = 0;
        tpl_ctrls->subpel_depth            = FULL_PEL;
        break;
    case 5:
        tpl_ctrls->enable                  = 1;
        tpl_ctrls->compute_rate            = 0;
        tpl_ctrls->enable_tpl_qps          = 0;
        tpl_ctrls->disable_intra_pred_nref = 1;
        tpl_ctrls->intra_mode_end          = DC_PRED;
        tpl_ctrls->reduced_tpl_group       = pcs->hierarchical_levels == 5
                  ? is_islice ? 4 : (resolution <= INPUT_SIZE_480p_RANGE ? 3 : 2)
                  : is_islice ? 3
                              : (resolution <= INPUT_SIZE_480p_RANGE ? 2 : 1);
        tpl_ctrls->pf_shape                = resolution <= INPUT_SIZE_480p_RANGE ? N2_SHAPE : N4_SHAPE;
        tpl_ctrls->use_sad_in_src_search   = 1;
        tpl_ctrls->dispenser_search_level  = 0;
        tpl_ctrls->subsample_tx            = 1;
        tpl_ctrls->subpel_depth            = FULL_PEL;
        break;
    case 6:
        tpl_ctrls->enable                  = 1;
        tpl_ctrls->compute_rate            = 0;
        tpl_ctrls->enable_tpl_qps          = 0;
        tpl_ctrls->disable_intra_pred_nref = 1;
        tpl_ctrls->intra_mode_end          = DC_PRED;
        tpl_ctrls->reduced_tpl_group       = pcs->hierarchical_levels == 5
                  ? is_islice ? 4 : (resolution <= INPUT_SIZE_480p_RANGE ? 3 : 2)
                  : pcs->hierarchical_levels == 4 ? is_islice ? 3 : (resolution <= INPUT_SIZE_480p_RANGE ? 2 : 1)
                  : is_islice                     ? 3
                                                  : (resolution <= INPUT_SIZE_480p_RANGE ? 2 : 0);
        tpl_ctrls->pf_shape                = resolution <= INPUT_SIZE_480p_RANGE ? N2_SHAPE : N4_SHAPE;
        tpl_ctrls->use_sad_in_src_search   = 1;
        tpl_ctrls->dispenser_search_level  = 1;
        tpl_ctrls->subsample_tx            = 2;
        tpl_ctrls->subpel_depth            = FULL_PEL;
        break;
    case 7:
        tpl_ctrls->enable                  = 1;
        tpl_ctrls->compute_rate            = 0;
        tpl_ctrls->enable_tpl_qps          = 0;
        tpl_ctrls->disable_intra_pred_nref = 1;
        tpl_ctrls->intra_mode_end          = DC_PRED;
        tpl_ctrls->reduced_tpl_group       = pcs->hierarchical_levels == 5
                  ? is_islice ? 4 : (resolution <= INPUT_SIZE_480p_RANGE ? 3 : 1)
                  : pcs->hierarchical_levels == 4 ? is_islice ? 3 : (resolution <= INPUT_SIZE_480p_RANGE ? 2 : 1)
                  : is_islice                     ? 3
                                                  : (resolution <= INPUT_SIZE_480p_RANGE ? 2 : 0);
        tpl_ctrls->pf_shape                = resolution <= INPUT_SIZE_480p_RANGE ? N2_SHAPE : N4_SHAPE;
        tpl_ctrls->use_sad_in_src_search   = 1;
        tpl_ctrls->dispenser_search_level  = 1;
        tpl_ctrls->subsample_tx            = 2;
        tpl_ctrls->subpel_depth            = FULL_PEL;
        break;
    case 8:
        tpl_ctrls->enable                  = 1;
        tpl_ctrls->compute_rate            = 0;
        tpl_ctrls->enable_tpl_qps          = 0;
        tpl_ctrls->disable_intra_pred_nref = 1;
        tpl_ctrls->intra_mode_end          = DC_PRED;
        tpl_ctrls->reduced_tpl_group       = pcs->hierarchical_levels == 5
                  ? is_islice ? 4 : (resolution <= INPUT_SIZE_480p_RANGE ? 3 : 1)
                  : is_islice ? 3
                              : (resolution <= INPUT_SIZE_480p_RANGE ? 2 : 0);
        tpl_ctrls->pf_shape                = resolution <= INPUT_SIZE_480p_RANGE ? N2_SHAPE : N4_SHAPE;
        tpl_ctrls->use_sad_in_src_search   = 1;
        tpl_ctrls->dispenser_search_level  = 1;
        tpl_ctrls->subsample_tx            = 2;
        tpl_ctrls->subpel_depth            = FULL_PEL;
        break;
    default: assert(0); break;
    }

    // Check user-defined settings for MAX intra mode
    if (scs->enable_paeth == 0)
        tpl_ctrls->intra_mode_end = MIN(tpl_ctrls->intra_mode_end, SMOOTH_H_PRED);

    if (scs->enable_smooth == 0)
        tpl_ctrls->intra_mode_end = MIN(tpl_ctrls->intra_mode_end, D67_PRED);

    // Derive synthesizer block size from frame size and tpl level
    tpl_ctrls->synth_blk_size = svt_aom_get_tpl_synthesizer_block_size(
        tpl_level, pcs->aligned_width, pcs->aligned_height);

    if ((int)pcs->hierarchical_levels <= tpl_ctrls->reduced_tpl_group)
        tpl_ctrls->reduced_tpl_group = -1;

    // TPL may only look at a subset of available pictures in tpl group, which may affect the r0 calcuation.
    // As a result, we defined a factor to adjust r0 (to compensate for TPL not using all available frames).
    if (tpl_ctrls->reduced_tpl_group >= 0) {
        switch ((pcs->hierarchical_levels - tpl_ctrls->reduced_tpl_group)) {
        case 0:
        default: tpl_ctrls->r0_adjust_factor = 0; break;
        case 1:
            tpl_ctrls->r0_adjust_factor = pcs->hierarchical_levels <= 2 ? 0.4
                : pcs->hierarchical_levels <= 3                         ? 0.8
                                                                        : 1.6;
            break;
        case 2:
            tpl_ctrls->r0_adjust_factor = pcs->hierarchical_levels <= 2 ? 0.6
                : pcs->hierarchical_levels <= 3                         ? 1.2
                                                                        : 2.4;
            break;
        case 3: tpl_ctrls->r0_adjust_factor = pcs->hierarchical_levels <= 3 ? 1.4 : 2.8; break;
        case 4: tpl_ctrls->r0_adjust_factor = 4.0; break;
        case 5: tpl_ctrls->r0_adjust_factor = 6.0; break;
        }

        // Adjust r0 scaling factor based on GOP structure and lookahead
        if (!scs->tpl_lad_mg)
            tpl_ctrls->r0_adjust_factor *= 1.5;
    } else {
        // No r0 adjustment when all frames are used
        tpl_ctrls->r0_adjust_factor = 0;

        // If no lookahead, apply r0 scaling
        if (!scs->tpl_lad_mg) {
            tpl_ctrls->r0_adjust_factor = is_islice ? 0
                : pcs->hierarchical_levels <= 2     ? 0.4
                : pcs->hierarchical_levels <= 3     ? 0.8
                                                    : 1.6;
        }
    }
}

static void set_cdef_controls(PictureParentControlSet *pcs, uint8_t cdef_level, Bool fast_decode) {
    CdefControls *cdef_ctrls = &pcs->cdef_ctrls;
    int           i, j, sf_idx, second_pass_fs_num;
    cdef_ctrls->use_reference_cdef_fs = 0;
    cdef_ctrls->use_skip_detector     = 0;
    switch (cdef_level) {
        // OFF
    case 0: cdef_ctrls->enabled = 0; break;
    case 1:
        // pf_set {0,1,..,15}
        // sf_set {0,1,2,3}
        cdef_ctrls->enabled                    = 1;
        cdef_ctrls->first_pass_fs_num          = 16;
        second_pass_fs_num                     = 3;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0]   = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1]   = pf_gi[1];
        cdef_ctrls->default_first_pass_fs[2]   = pf_gi[2];
        cdef_ctrls->default_first_pass_fs[3]   = pf_gi[3];
        cdef_ctrls->default_first_pass_fs[4]   = pf_gi[4];
        cdef_ctrls->default_first_pass_fs[5]   = pf_gi[5];
        cdef_ctrls->default_first_pass_fs[6]   = pf_gi[6];
        cdef_ctrls->default_first_pass_fs[7]   = pf_gi[7];
        cdef_ctrls->default_first_pass_fs[8]   = pf_gi[8];
        cdef_ctrls->default_first_pass_fs[9]   = pf_gi[9];
        cdef_ctrls->default_first_pass_fs[10]  = pf_gi[10];
        cdef_ctrls->default_first_pass_fs[11]  = pf_gi[11];
        cdef_ctrls->default_first_pass_fs[12]  = pf_gi[12];
        cdef_ctrls->default_first_pass_fs[13]  = pf_gi[13];
        cdef_ctrls->default_first_pass_fs[14]  = pf_gi[14];
        cdef_ctrls->default_first_pass_fs[15]  = pf_gi[15];
        sf_idx                                 = 0;
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++) {
            int pf_idx = cdef_ctrls->default_first_pass_fs[i];
            for (j = 1; j < 4; j++) {
                cdef_ctrls->default_second_pass_fs[sf_idx] = pf_idx + j;
                sf_idx++;
            }
        }
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++)
            cdef_ctrls->default_first_pass_fs_uv[i] = cdef_ctrls->default_first_pass_fs[i];
        for (i = 0; i < cdef_ctrls->default_second_pass_fs_num; i++)
            cdef_ctrls->default_second_pass_fs_uv[i] = cdef_ctrls->default_second_pass_fs[i];
        cdef_ctrls->use_reference_cdef_fs = 0;
        cdef_ctrls->search_best_ref_fs    = 0;
        cdef_ctrls->subsampling_factor    = 1;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 2: // N
        // pf_set {0,1,2,4,5,6,8,9,10,12,13,14}
        // sf_set {0,1,..,3}
        cdef_ctrls->enabled                    = 1;
        cdef_ctrls->first_pass_fs_num          = 12;
        second_pass_fs_num                     = 3;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0]   = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1]   = pf_gi[1];
        cdef_ctrls->default_first_pass_fs[2]   = pf_gi[2];
        cdef_ctrls->default_first_pass_fs[3]   = pf_gi[4];
        cdef_ctrls->default_first_pass_fs[4]   = pf_gi[5];
        cdef_ctrls->default_first_pass_fs[5]   = pf_gi[6];
        cdef_ctrls->default_first_pass_fs[6]   = pf_gi[8];
        cdef_ctrls->default_first_pass_fs[7]   = pf_gi[9];
        cdef_ctrls->default_first_pass_fs[8]   = pf_gi[10];
        cdef_ctrls->default_first_pass_fs[9]   = pf_gi[12];
        cdef_ctrls->default_first_pass_fs[10]  = pf_gi[13];
        cdef_ctrls->default_first_pass_fs[11]  = pf_gi[14];
        sf_idx                                 = 0;
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++) {
            int pf_idx = cdef_ctrls->default_first_pass_fs[i];
            for (j = 1; j < 4; j++) {
                cdef_ctrls->default_second_pass_fs[sf_idx] = pf_idx + j;
                sf_idx++;
            }
        }
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++)
            cdef_ctrls->default_first_pass_fs_uv[i] = cdef_ctrls->default_first_pass_fs[i];
        for (i = 0; i < cdef_ctrls->default_second_pass_fs_num; i++)
            cdef_ctrls->default_second_pass_fs_uv[i] = -1; // cdef_ctrls->default_second_pass_fs[i];
        cdef_ctrls->use_reference_cdef_fs = 0;
        cdef_ctrls->search_best_ref_fs    = 0;
        cdef_ctrls->subsampling_factor    = 1;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 3:
        // pf_set {0,2,4,6,8,10,12,14}
        // sf_set {0,1,..,3}
        cdef_ctrls->enabled                    = 1;
        cdef_ctrls->first_pass_fs_num          = 8;
        second_pass_fs_num                     = 3;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0]   = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1]   = pf_gi[2];
        cdef_ctrls->default_first_pass_fs[2]   = pf_gi[4];
        cdef_ctrls->default_first_pass_fs[3]   = pf_gi[6];
        cdef_ctrls->default_first_pass_fs[4]   = pf_gi[8];
        cdef_ctrls->default_first_pass_fs[5]   = pf_gi[10];
        cdef_ctrls->default_first_pass_fs[6]   = pf_gi[12];
        cdef_ctrls->default_first_pass_fs[7]   = pf_gi[14];
        sf_idx                                 = 0;
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++) {
            int pf_idx = cdef_ctrls->default_first_pass_fs[i];
            for (j = 1; j < 4; j++) {
                cdef_ctrls->default_second_pass_fs[sf_idx] = pf_idx + j;
                sf_idx++;
            }
        }
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++)
            cdef_ctrls->default_first_pass_fs_uv[i] = cdef_ctrls->default_first_pass_fs[i];
        for (i = 0; i < cdef_ctrls->default_second_pass_fs_num; i++)
            cdef_ctrls->default_second_pass_fs_uv[i] = cdef_ctrls->default_second_pass_fs[i];
        cdef_ctrls->use_reference_cdef_fs = 0;
        cdef_ctrls->search_best_ref_fs    = 0;
        cdef_ctrls->subsampling_factor    = 1;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 4:
        // pf_set {0,4,8,12,15}
        // sf_set {0,1,..,3}
        cdef_ctrls->enabled                    = 1;
        cdef_ctrls->first_pass_fs_num          = 5;
        second_pass_fs_num                     = 3;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0]   = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1]   = pf_gi[4];
        cdef_ctrls->default_first_pass_fs[2]   = pf_gi[8];
        cdef_ctrls->default_first_pass_fs[3]   = pf_gi[12];
        cdef_ctrls->default_first_pass_fs[4]   = pf_gi[15];
        sf_idx                                 = 0;
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++) {
            int pf_idx = cdef_ctrls->default_first_pass_fs[i];
            for (j = 1; j < 4; j++) {
                cdef_ctrls->default_second_pass_fs[sf_idx] = pf_idx + j;
                sf_idx++;
            }
        }
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++)
            cdef_ctrls->default_first_pass_fs_uv[i] = cdef_ctrls->default_first_pass_fs[i];
        for (i = 0; i < cdef_ctrls->default_second_pass_fs_num; i++)
            cdef_ctrls->default_second_pass_fs_uv[i] = -1; // cdef_ctrls->default_second_pass_fs[i];
        cdef_ctrls->use_reference_cdef_fs = 0;
        cdef_ctrls->search_best_ref_fs    = 0;
        cdef_ctrls->subsampling_factor    = 1;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 5:
        // pf_set {0,5,10,15}
        // sf_set {0,1,..,3}
        cdef_ctrls->enabled                    = 1;
        cdef_ctrls->first_pass_fs_num          = 4;
        second_pass_fs_num                     = 3;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0]   = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1]   = pf_gi[5];
        cdef_ctrls->default_first_pass_fs[2]   = pf_gi[10];
        cdef_ctrls->default_first_pass_fs[3]   = pf_gi[15];
        sf_idx                                 = 0;
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++) {
            int pf_idx = cdef_ctrls->default_first_pass_fs[i];
            for (j = 1; j < 4; j++) {
                cdef_ctrls->default_second_pass_fs[sf_idx] = pf_idx + j;
                sf_idx++;
            }
        }
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++)
            cdef_ctrls->default_first_pass_fs_uv[i] = cdef_ctrls->default_first_pass_fs[i];
        for (i = 0; i < cdef_ctrls->default_second_pass_fs_num; i++)
            cdef_ctrls->default_second_pass_fs_uv[i] = cdef_ctrls->default_second_pass_fs[i];
        cdef_ctrls->use_reference_cdef_fs = 0;
        cdef_ctrls->search_best_ref_fs    = 0;
        cdef_ctrls->subsampling_factor    = 1;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 6:
        // pf_set {0,7,15}
        // sf_set {0,1,..,3}
        cdef_ctrls->enabled                    = 1;
        cdef_ctrls->first_pass_fs_num          = 3;
        second_pass_fs_num                     = 3;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0]   = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1]   = pf_gi[7];
        cdef_ctrls->default_first_pass_fs[2]   = pf_gi[15];
        sf_idx                                 = 0;
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++) {
            int pf_idx = cdef_ctrls->default_first_pass_fs[i];
            for (j = 1; j < 4; j++) {
                cdef_ctrls->default_second_pass_fs[sf_idx] = pf_idx + j;
                sf_idx++;
            }
        }
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++)
            cdef_ctrls->default_first_pass_fs_uv[i] = cdef_ctrls->default_first_pass_fs[i];
        for (i = 0; i < cdef_ctrls->default_second_pass_fs_num; i++)
            cdef_ctrls->default_second_pass_fs_uv[i] = -1; // cdef_ctrls->default_second_pass_fs[i];
        cdef_ctrls->use_reference_cdef_fs = 0;
        cdef_ctrls->search_best_ref_fs    = 0;
        cdef_ctrls->subsampling_factor    = 1;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 7:
        // pf_set {0,7,15}
        // sf_set {0,1,2}
        cdef_ctrls->enabled                    = 1;
        cdef_ctrls->first_pass_fs_num          = 3;
        second_pass_fs_num                     = 2;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0]   = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1]   = pf_gi[7];
        cdef_ctrls->default_first_pass_fs[2]   = pf_gi[15];

        cdef_ctrls->default_second_pass_fs[0] = pf_gi[0] + 1;
        cdef_ctrls->default_second_pass_fs[1] = pf_gi[0] + 2;
        cdef_ctrls->default_second_pass_fs[2] = pf_gi[7] + 1;
        cdef_ctrls->default_second_pass_fs[3] = pf_gi[7] + 2;
        cdef_ctrls->default_second_pass_fs[4] = pf_gi[15] + 1;
        cdef_ctrls->default_second_pass_fs[5] = pf_gi[15] + 2;
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++)
            cdef_ctrls->default_first_pass_fs_uv[i] = cdef_ctrls->default_first_pass_fs[i];
        for (i = 0; i < cdef_ctrls->default_second_pass_fs_num; i++)
            cdef_ctrls->default_second_pass_fs_uv[i] = cdef_ctrls->default_second_pass_fs[i];
        cdef_ctrls->use_reference_cdef_fs = 0;
        cdef_ctrls->search_best_ref_fs    = 0;
        cdef_ctrls->subsampling_factor    = 1;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 8:
        // pf_set {0,7,15}
        // sf_set {0,2}
        cdef_ctrls->enabled                    = 1;
        cdef_ctrls->first_pass_fs_num          = 3;
        second_pass_fs_num                     = 1;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0]   = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1]   = pf_gi[7];
        cdef_ctrls->default_first_pass_fs[2]   = pf_gi[15];

        cdef_ctrls->default_second_pass_fs[0] = pf_gi[0] + 2;
        cdef_ctrls->default_second_pass_fs[1] = pf_gi[7] + 2;
        cdef_ctrls->default_second_pass_fs[2] = pf_gi[15] + 2;
        for (i = 0; i < cdef_ctrls->first_pass_fs_num; i++)
            cdef_ctrls->default_first_pass_fs_uv[i] = cdef_ctrls->default_first_pass_fs[i];
        for (i = 0; i < cdef_ctrls->default_second_pass_fs_num; i++)
            cdef_ctrls->default_second_pass_fs_uv[i] = -1; // cdef_ctrls->default_second_pass_fs[i];
        cdef_ctrls->use_reference_cdef_fs = 0;
        cdef_ctrls->search_best_ref_fs    = 0;
        cdef_ctrls->subsampling_factor    = 1;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 9:
        // pf_set {0,15}
        // sf_set {0,2}
        cdef_ctrls->enabled                    = 1;
        cdef_ctrls->first_pass_fs_num          = 2;
        second_pass_fs_num                     = 1;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0]   = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1]   = pf_gi[15];

        cdef_ctrls->default_second_pass_fs[0]    = pf_gi[0] + 2;
        cdef_ctrls->default_second_pass_fs[1]    = pf_gi[15] + 2;
        cdef_ctrls->default_first_pass_fs_uv[0]  = cdef_ctrls->default_first_pass_fs[0];
        cdef_ctrls->default_first_pass_fs_uv[1]  = cdef_ctrls->default_first_pass_fs[1];
        cdef_ctrls->default_first_pass_fs_uv[2]  = -1; // if using search_best_ref_fs, set at least 3 filters
        cdef_ctrls->default_second_pass_fs_uv[0] = -1; // cdef_ctrls->default_second_pass_fs[0];
        cdef_ctrls->default_second_pass_fs_uv[1] = -1; // cdef_ctrls->default_second_pass_fs[1];
        cdef_ctrls->use_reference_cdef_fs        = 0;
        cdef_ctrls->search_best_ref_fs           = 0;
        cdef_ctrls->subsampling_factor           = 4;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 10:
        // pf_set {0,15}
        // sf_set {0,2}
        cdef_ctrls->enabled                    = 1;
        cdef_ctrls->first_pass_fs_num          = 2;
        second_pass_fs_num                     = 1;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0]   = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1]   = pf_gi[15];

        cdef_ctrls->default_second_pass_fs[0] = pf_gi[0] + 2;
        cdef_ctrls->default_second_pass_fs[1] = pf_gi[15] + 2;

        cdef_ctrls->default_first_pass_fs_uv[0]  = cdef_ctrls->default_first_pass_fs[0];
        cdef_ctrls->default_first_pass_fs_uv[1]  = cdef_ctrls->default_first_pass_fs[1];
        cdef_ctrls->default_first_pass_fs_uv[2]  = -1; // when using search_best_ref_fs, set at least 3 filters
        cdef_ctrls->default_second_pass_fs_uv[0] = -1;
        cdef_ctrls->default_second_pass_fs_uv[1] = -1;

        cdef_ctrls->use_reference_cdef_fs = 0;
        cdef_ctrls->search_best_ref_fs    = 1;
        cdef_ctrls->subsampling_factor    = 4;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 11:
        // pf_set {0,15}
        // sf_set {0,2}
        cdef_ctrls->enabled                    = 1;
        cdef_ctrls->first_pass_fs_num          = 2;
        second_pass_fs_num                     = 1;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0]   = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1]   = pf_gi[15];

        cdef_ctrls->default_second_pass_fs[0] = pf_gi[0] + 2;
        cdef_ctrls->default_second_pass_fs[1] = pf_gi[15] + 2;

        cdef_ctrls->default_first_pass_fs_uv[0]  = cdef_ctrls->default_first_pass_fs[0];
        cdef_ctrls->default_first_pass_fs_uv[1]  = cdef_ctrls->default_first_pass_fs[1];
        cdef_ctrls->default_first_pass_fs_uv[2]  = -1; // if using search_best_ref_fs, set at least 3 filters
        cdef_ctrls->default_second_pass_fs_uv[0] = -1;
        cdef_ctrls->default_second_pass_fs_uv[1] = -1;
        cdef_ctrls->use_reference_cdef_fs        = 1;
        cdef_ctrls->search_best_ref_fs           = 1;
        cdef_ctrls->subsampling_factor           = 4;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 12:
        // pf_set {0,15}
        // sf_set {0,2}
        cdef_ctrls->enabled                    = 1;
        cdef_ctrls->first_pass_fs_num          = 2;
        second_pass_fs_num                     = 1;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0]   = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1]   = pf_gi[15];

        cdef_ctrls->default_second_pass_fs[0]    = pf_gi[0] + 2;
        cdef_ctrls->default_second_pass_fs[1]    = pf_gi[15] + 2;
        cdef_ctrls->default_first_pass_fs_uv[0]  = cdef_ctrls->default_first_pass_fs[0];
        cdef_ctrls->default_first_pass_fs_uv[1]  = cdef_ctrls->default_first_pass_fs[1];
        cdef_ctrls->default_first_pass_fs_uv[2]  = -1; // if using search_best_ref_fs, set at least 3 filters
        cdef_ctrls->default_second_pass_fs_uv[0] = -1; // cdef_ctrls->default_second_pass_fs[0];
        cdef_ctrls->default_second_pass_fs_uv[1] = -1; // cdef_ctrls->default_second_pass_fs[1];
        cdef_ctrls->use_reference_cdef_fs        = 0;
        cdef_ctrls->search_best_ref_fs           = 0;
        cdef_ctrls->subsampling_factor           = 4;
        cdef_ctrls->use_skip_detector            = 0;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 13:
        // pf_set {0,15}
        // sf_set {0,2}
        cdef_ctrls->enabled                    = 1;
        cdef_ctrls->first_pass_fs_num          = 2;
        second_pass_fs_num                     = 1;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0]   = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1]   = pf_gi[15];

        cdef_ctrls->default_second_pass_fs[0] = pf_gi[0] + 2;
        cdef_ctrls->default_second_pass_fs[1] = pf_gi[15] + 2;

        cdef_ctrls->default_first_pass_fs_uv[0]  = cdef_ctrls->default_first_pass_fs[0];
        cdef_ctrls->default_first_pass_fs_uv[1]  = cdef_ctrls->default_first_pass_fs[1];
        cdef_ctrls->default_first_pass_fs_uv[2]  = -1; // when using search_best_ref_fs, set at least 3 filters
        cdef_ctrls->default_second_pass_fs_uv[0] = -1;
        cdef_ctrls->default_second_pass_fs_uv[1] = -1;

        cdef_ctrls->use_reference_cdef_fs = 0;
        cdef_ctrls->search_best_ref_fs    = 1;
        cdef_ctrls->subsampling_factor    = 4;
        cdef_ctrls->use_skip_detector     = 0;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 14:
        // pf_set {0,15}
        // sf_set {0,2}
        cdef_ctrls->enabled                    = 1;
        cdef_ctrls->first_pass_fs_num          = 2;
        second_pass_fs_num                     = 1;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0]   = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1]   = pf_gi[15];

        cdef_ctrls->default_second_pass_fs[0] = pf_gi[0] + 2;
        cdef_ctrls->default_second_pass_fs[1] = pf_gi[15] + 2;

        cdef_ctrls->default_first_pass_fs_uv[0]  = cdef_ctrls->default_first_pass_fs[0];
        cdef_ctrls->default_first_pass_fs_uv[1]  = cdef_ctrls->default_first_pass_fs[1];
        cdef_ctrls->default_first_pass_fs_uv[2]  = -1; // if using search_best_ref_fs, set at least 3 filters
        cdef_ctrls->default_second_pass_fs_uv[0] = -1;
        cdef_ctrls->default_second_pass_fs_uv[1] = -1;
        cdef_ctrls->use_reference_cdef_fs        = 1;
        cdef_ctrls->search_best_ref_fs           = 1;
        cdef_ctrls->subsampling_factor           = 4;
        cdef_ctrls->use_skip_detector            = 0;
        if (fast_decode == 0)
            cdef_ctrls->zero_fs_cost_bias = 0;
        else
            cdef_ctrls->zero_fs_cost_bias = pcs->input_resolution <= INPUT_SIZE_360p_RANGE ? 0 : 62;
        break;
    case 15:
        // pf_set {0,15}
        // sf_set {0,2}
        cdef_ctrls->enabled                    = 1;
        cdef_ctrls->first_pass_fs_num          = 2;
        second_pass_fs_num                     = 1;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0]   = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1]   = pf_gi[15];

        cdef_ctrls->default_second_pass_fs[0]    = pf_gi[0] + 2;
        cdef_ctrls->default_second_pass_fs[1]    = pf_gi[15] + 2;
        cdef_ctrls->default_first_pass_fs_uv[0]  = cdef_ctrls->default_first_pass_fs[0];
        cdef_ctrls->default_first_pass_fs_uv[1]  = cdef_ctrls->default_first_pass_fs[1];
        cdef_ctrls->default_first_pass_fs_uv[2]  = -1; // if using search_best_ref_fs, set at least 3 filters
        cdef_ctrls->default_second_pass_fs_uv[0] = -1; // cdef_ctrls->default_second_pass_fs[0];
        cdef_ctrls->default_second_pass_fs_uv[1] = -1; // cdef_ctrls->default_second_pass_fs[1];
        cdef_ctrls->use_reference_cdef_fs        = 0;
        cdef_ctrls->search_best_ref_fs           = 0;
        cdef_ctrls->subsampling_factor           = 4;
        cdef_ctrls->zero_fs_cost_bias            = 62;
        cdef_ctrls->use_skip_detector            = 0;
        break;
    case 16:
        // pf_set {0,15}
        // sf_set {0,2}
        cdef_ctrls->enabled                    = 1;
        cdef_ctrls->first_pass_fs_num          = 2;
        second_pass_fs_num                     = 1;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0]   = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1]   = pf_gi[15];

        cdef_ctrls->default_second_pass_fs[0] = pf_gi[0] + 2;
        cdef_ctrls->default_second_pass_fs[1] = pf_gi[15] + 2;

        cdef_ctrls->default_first_pass_fs_uv[0]  = cdef_ctrls->default_first_pass_fs[0];
        cdef_ctrls->default_first_pass_fs_uv[1]  = cdef_ctrls->default_first_pass_fs[1];
        cdef_ctrls->default_first_pass_fs_uv[2]  = -1; // when using search_best_ref_fs, set at least 3 filters
        cdef_ctrls->default_second_pass_fs_uv[0] = -1;
        cdef_ctrls->default_second_pass_fs_uv[1] = -1;

        cdef_ctrls->use_reference_cdef_fs = 0;
        cdef_ctrls->search_best_ref_fs    = 1;
        cdef_ctrls->subsampling_factor    = 4;
        cdef_ctrls->zero_fs_cost_bias     = 62;
        cdef_ctrls->use_skip_detector     = 1;
        break;
    case 17:
        // pf_set {0,15}
        // sf_set {0,2}
        cdef_ctrls->enabled                    = 1;
        cdef_ctrls->first_pass_fs_num          = 2;
        second_pass_fs_num                     = 1;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0]   = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1]   = pf_gi[15];

        cdef_ctrls->default_second_pass_fs[0] = pf_gi[0] + 2;
        cdef_ctrls->default_second_pass_fs[1] = pf_gi[15] + 2;

        cdef_ctrls->default_first_pass_fs_uv[0]  = cdef_ctrls->default_first_pass_fs[0];
        cdef_ctrls->default_first_pass_fs_uv[1]  = cdef_ctrls->default_first_pass_fs[1];
        cdef_ctrls->default_first_pass_fs_uv[2]  = -1; // if using search_best_ref_fs, set at least 3 filters
        cdef_ctrls->default_second_pass_fs_uv[0] = -1;
        cdef_ctrls->default_second_pass_fs_uv[1] = -1;
        cdef_ctrls->use_reference_cdef_fs        = 1;
        cdef_ctrls->search_best_ref_fs           = 1;
        cdef_ctrls->subsampling_factor           = 4;
        cdef_ctrls->zero_fs_cost_bias            = 62;
        cdef_ctrls->use_skip_detector            = 1;
        break;
    case 18:
        // pf_set {0}
        // sf_set {0}
        cdef_ctrls->enabled                     = 1;
        cdef_ctrls->first_pass_fs_num           = 1;
        second_pass_fs_num                      = 0;
        cdef_ctrls->default_second_pass_fs_num  = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0]    = pf_gi[0];
        cdef_ctrls->default_first_pass_fs_uv[0] = cdef_ctrls->default_first_pass_fs[0];
        cdef_ctrls->use_reference_cdef_fs       = 1;
        cdef_ctrls->search_best_ref_fs          = 1;
        cdef_ctrls->subsampling_factor          = 4;
        cdef_ctrls->zero_fs_cost_bias           = 62;
        cdef_ctrls->use_skip_detector           = 1;
        break;
    case 19:
        // pf_set {0,15}
        // sf_set {0,2}
        cdef_ctrls->enabled                    = 1;
        cdef_ctrls->first_pass_fs_num          = 2;
        second_pass_fs_num                     = 1;
        cdef_ctrls->default_second_pass_fs_num = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0]   = pf_gi[0];
        cdef_ctrls->default_first_pass_fs[1]   = pf_gi[15];

        cdef_ctrls->default_second_pass_fs[0]    = pf_gi[0] + 2;
        cdef_ctrls->default_second_pass_fs[1]    = pf_gi[15] + 2;
        cdef_ctrls->default_first_pass_fs_uv[0]  = cdef_ctrls->default_first_pass_fs[0];
        cdef_ctrls->default_first_pass_fs_uv[1]  = cdef_ctrls->default_first_pass_fs[1];
        cdef_ctrls->default_second_pass_fs_uv[0] = -1;
        cdef_ctrls->default_second_pass_fs_uv[1] = -1;
        cdef_ctrls->use_reference_cdef_fs        = 0;
        cdef_ctrls->search_best_ref_fs           = 0;
        cdef_ctrls->subsampling_factor           = 4;
        cdef_ctrls->zero_fs_cost_bias            = 0;
        cdef_ctrls->use_skip_detector            = 0;
        break;
    case 20:
        // pf_set {0}
        // sf_set {0}
        cdef_ctrls->enabled                     = 1;
        cdef_ctrls->first_pass_fs_num           = 1;
        second_pass_fs_num                      = 0;
        cdef_ctrls->default_second_pass_fs_num  = cdef_ctrls->first_pass_fs_num * second_pass_fs_num;
        cdef_ctrls->default_first_pass_fs[0]    = pf_gi[0];
        cdef_ctrls->default_first_pass_fs_uv[0] = cdef_ctrls->default_first_pass_fs[0];
        cdef_ctrls->use_reference_cdef_fs       = 1;
        cdef_ctrls->search_best_ref_fs          = 0;
        cdef_ctrls->subsampling_factor          = 4;
        cdef_ctrls->zero_fs_cost_bias           = 0;
        cdef_ctrls->use_skip_detector           = 0;
        break;
    default: assert(0); break;
    }
}

static void svt_aom_set_wn_filter_ctrls(Av1Common *cm, uint8_t wn_filter_lvl) {
    WnFilterCtrls *ctrls = &cm->wn_filter_ctrls;

    switch (wn_filter_lvl) {
    case 0: ctrls->enabled = 0; break;
    case 1:
        ctrls->enabled                 = 1;
        ctrls->use_chroma              = 1;
        ctrls->filter_tap_lvl          = 1;
        ctrls->use_refinement          = 1;
        ctrls->max_one_refinement_step = 0;
        ctrls->use_prev_frame_coeffs   = 0;
        break;
    case 2:
        ctrls->enabled                 = 1;
        ctrls->use_chroma              = 1;
        ctrls->filter_tap_lvl          = 1;
        ctrls->use_refinement          = 1;
        ctrls->max_one_refinement_step = 1;
        ctrls->use_prev_frame_coeffs   = 0;
        break;
    case 3:
        ctrls->enabled                 = 1;
        ctrls->use_chroma              = 1;
        ctrls->filter_tap_lvl          = 2;
        ctrls->use_refinement          = 1;
        ctrls->max_one_refinement_step = 1;
        ctrls->use_prev_frame_coeffs   = 0;
        break;
    case 4:
        ctrls->enabled                 = 1;
        ctrls->use_chroma              = 1;
        ctrls->filter_tap_lvl          = 2;
        ctrls->use_refinement          = 0;
        ctrls->max_one_refinement_step = 1;
        ctrls->use_prev_frame_coeffs   = 0;
        break;
    case 5:
        ctrls->enabled                 = 1;
        ctrls->use_chroma              = 0;
        ctrls->filter_tap_lvl          = 2;
        ctrls->use_refinement          = 0;
        ctrls->max_one_refinement_step = 1;
        ctrls->use_prev_frame_coeffs   = 0;
        break;
    case 6:
        ctrls->enabled                 = 1;
        ctrls->use_chroma              = 0;
        ctrls->filter_tap_lvl          = 2;
        ctrls->use_refinement          = 0;
        ctrls->max_one_refinement_step = 1;
        ctrls->use_prev_frame_coeffs   = 1;
        break;
    default: assert(0); break;
    }
}

static void svt_aom_set_sg_filter_ctrls(Av1Common *cm, uint8_t sg_filter_lvl) {
    SgFilterCtrls *ctrls = &cm->sg_filter_ctrls;

    switch (sg_filter_lvl) {
    case 0: ctrls->enabled = 0; break;
    case 1:
        ctrls->enabled     = 1;
        ctrls->use_chroma  = 1;
        ctrls->step_range  = 16;
        ctrls->start_ep[0] = 0;
        ctrls->end_ep[0]   = 16;
        ctrls->ep_inc[0]   = 1;
        ctrls->start_ep[1] = 0;
        ctrls->end_ep[1]   = 16;
        ctrls->ep_inc[1]   = 1;
        ctrls->refine[0]   = 1;
        ctrls->refine[1]   = 1;
        break;
    case 2:
        ctrls->enabled     = 1;
        ctrls->use_chroma  = 1;
        ctrls->step_range  = 16;
        ctrls->start_ep[0] = 0;
        ctrls->end_ep[0]   = 16;
        ctrls->ep_inc[0]   = 1;
        ctrls->start_ep[1] = 4;
        ctrls->end_ep[1]   = 5;
        ctrls->ep_inc[1]   = 1;
        ctrls->refine[0]   = 1;
        ctrls->refine[1]   = 0;
        break;
    case 3:
        ctrls->enabled     = 1;
        ctrls->use_chroma  = 1;
        ctrls->step_range  = 16;
        ctrls->start_ep[0] = 0;
        ctrls->end_ep[0]   = 16;
        ctrls->ep_inc[0]   = 8;
        ctrls->start_ep[1] = 4;
        ctrls->end_ep[1]   = 5;
        ctrls->ep_inc[1]   = 1;
        ctrls->refine[0]   = 1;
        ctrls->refine[1]   = 0;
        break;
    case 4:
        ctrls->enabled     = 1;
        ctrls->use_chroma  = 1;
        ctrls->step_range  = 0;
        ctrls->start_ep[0] = 0;
        ctrls->end_ep[0]   = 16;
        ctrls->ep_inc[0]   = 1;
        ctrls->start_ep[1] = 0;
        ctrls->end_ep[1]   = 16;
        ctrls->ep_inc[1]   = 1;
        ctrls->refine[0]   = 1;
        ctrls->refine[1]   = 1;
        break;
    default: assert(0); break;
    }
}

// Returns the level for Wiener restoration filter
static uint8_t svt_aom_get_wn_filter_level(EncMode enc_mode, uint8_t input_resolution, Bool is_ref) {
    uint8_t wn_filter_lvl = 0;
    if (enc_mode <= ENC_M2)
        wn_filter_lvl = 1;
    else if (enc_mode <= ENC_M9)
        wn_filter_lvl = is_ref ? 5 : 0;
    else
        wn_filter_lvl = 0;

    // higher resolutions will shut restoration to save memory
    if (input_resolution >= INPUT_SIZE_8K_RANGE)
        wn_filter_lvl = 0;

    return wn_filter_lvl;
}

// Returns the level for self-guided restoration filter
static uint8_t svt_aom_get_sg_filter_level(EncMode enc_mode, Bool fast_decode, uint8_t input_resolution, Bool is_base) {
    uint8_t sg_filter_lvl = 0;
    if (fast_decode == 0) {
        if (enc_mode <= ENC_M0)
            sg_filter_lvl = 1;
        else if (enc_mode <= ENC_M4)
            sg_filter_lvl = 3;
        else
            sg_filter_lvl = 0;
    } else {
        if (enc_mode <= ENC_M1)
            sg_filter_lvl = input_resolution <= INPUT_SIZE_360p_RANGE ? 1 : 0;
        else if (enc_mode <= ENC_M4)
            sg_filter_lvl = input_resolution <= INPUT_SIZE_360p_RANGE ? (is_base ? 1 : 3) : 0;
        else if (enc_mode <= ENC_M5)
            sg_filter_lvl = input_resolution <= INPUT_SIZE_360p_RANGE ? (is_base ? 1 : 0) : 0;
        else
            sg_filter_lvl = 0;
    }

    // higher resolutions will shut restoration to save memory
    if (input_resolution >= INPUT_SIZE_8K_RANGE)
        sg_filter_lvl = 0;

    return sg_filter_lvl;
}
/*
* return the DLF level
* Used by svt_aom_sig_deriv_multi_processes and memory allocation
*
* 0 -- DLF OFF
* 1 -- Full Frame DLF
* 2 -- SB-Based DLF
* 3 -- SB-Based DLF + skip DLF if SB has 0 coeffs
*/
static uint8_t get_dlf_level(EncMode enc_mode, uint8_t is_ref, uint8_t is_16bit, Bool fast_decode,
                             EbInputResolution resolution, bool rtc_tune, uint8_t sc_class1) {
    uint8_t dlf_level;
    if (rtc_tune) {
        if (sc_class1) {
            if (enc_mode <= ENC_M10)
                dlf_level = resolution <= INPUT_SIZE_360p_RANGE ? 3 : (is_ref ? 4 : 0);
            else
                dlf_level = is_ref ? 5 : 0;
        } else if (enc_mode <= ENC_M9)
            dlf_level = resolution <= INPUT_SIZE_360p_RANGE ? 3 : (is_ref ? 3 : 5);
        else if (enc_mode <= ENC_M11)
            dlf_level = resolution <= INPUT_SIZE_360p_RANGE ? 3 : (is_ref ? 4 : 0);
        else
            dlf_level = is_ref ? 5 : 0;
    }
    // Don't disable DLF for low resolutions when fast-decode is used
    else if (fast_decode == 0 || resolution <= INPUT_SIZE_360p_RANGE) {
        if (enc_mode <= ENC_M2)
            dlf_level = 1;
        else if (enc_mode <= ENC_M6)
            dlf_level = 2;
        else if (enc_mode <= ENC_M7) {
            if (resolution <= INPUT_SIZE_1080p_RANGE)
                dlf_level = 3;
            else
                dlf_level = 2;
        } else if (enc_mode <= ENC_M8)
            dlf_level = resolution <= INPUT_SIZE_360p_RANGE ? 3 : (is_ref ? 3 : 5);
        else if (enc_mode <= ENC_M11)
            dlf_level = resolution <= INPUT_SIZE_360p_RANGE ? 3 : (is_ref ? 4 : 0);
        else if (enc_mode <= ENC_M12)
            dlf_level = is_ref ? 5 : 0;
        else
            dlf_level = (is_16bit && is_ref) ? 5 : 0;
    } else {
        dlf_level = is_ref ? 5 : 0;
    }
    return dlf_level;
}

static void svt_aom_set_dlf_controls(PictureParentControlSet *pcs, uint8_t dlf_level) {
    DlfCtrls *ctrls = &pcs->dlf_ctrls;

    switch (dlf_level) {
    case 0:
        ctrls->enabled                  = 0;
        ctrls->sb_based_dlf             = 0;
        ctrls->dlf_avg                  = 0;
        ctrls->dlf_avg_uv               = 0;
        ctrls->early_exit_convergence   = 0;
        ctrls->zero_filter_strength_lvl = 0;
        break;
    case 1:
        ctrls->enabled                  = 1;
        ctrls->sb_based_dlf             = 0;
        ctrls->dlf_avg                  = 0;
        ctrls->dlf_avg_uv               = 0;
        ctrls->early_exit_convergence   = 0;
        ctrls->zero_filter_strength_lvl = 0;
        break;
    case 2:
        ctrls->enabled                  = 1;
        ctrls->sb_based_dlf             = 0;
        ctrls->dlf_avg                  = 1;
        ctrls->dlf_avg_uv               = 1;
        ctrls->early_exit_convergence   = 1;
        ctrls->zero_filter_strength_lvl = 0;
        break;
    case 3:
        ctrls->enabled      = 1;
        ctrls->sb_based_dlf = 1;

        ctrls->dlf_avg                  = 0;
        ctrls->dlf_avg_uv               = 0;
        ctrls->early_exit_convergence   = 0;
        ctrls->zero_filter_strength_lvl = 1;
        break;
    case 4:
        ctrls->enabled                  = 1;
        ctrls->sb_based_dlf             = 1;
        ctrls->dlf_avg                  = 0;
        ctrls->dlf_avg_uv               = 0;
        ctrls->early_exit_convergence   = 0;
        ctrls->zero_filter_strength_lvl = 2;
        break;
    case 5:
        ctrls->enabled                  = 1;
        ctrls->sb_based_dlf             = 1;
        ctrls->dlf_avg                  = 0;
        ctrls->dlf_avg_uv               = 0;
        ctrls->early_exit_convergence   = 0;
        ctrls->zero_filter_strength_lvl = 3;
        break;
    default: assert(0); break;
    }
}

/*
    set controls for intra block copy
*/
static void set_intrabc_level(PictureParentControlSet *pcs, SequenceControlSet *scs, uint8_t ibc_level) {
    IntraBCCtrls *intraBC_ctrls = &pcs->intraBC_ctrls;

    switch (ibc_level) {
    case 0: intraBC_ctrls->enabled = 0; break;
    case 1:
        intraBC_ctrls->enabled             = pcs->sc_class1;
        intraBC_ctrls->ibc_shift           = 0;
        intraBC_ctrls->ibc_direction       = 0;
        intraBC_ctrls->hash_4x4_blocks     = !svt_aom_get_disallow_4x4(pcs->enc_mode, pcs->slice_type);
        intraBC_ctrls->max_block_size_hash = scs->super_block_size;
        break;
    case 2:
        intraBC_ctrls->enabled             = pcs->sc_class1;
        intraBC_ctrls->ibc_shift           = 1;
        intraBC_ctrls->ibc_direction       = 0;
        intraBC_ctrls->hash_4x4_blocks     = !svt_aom_get_disallow_4x4(pcs->enc_mode, pcs->slice_type);
        intraBC_ctrls->max_block_size_hash = scs->super_block_size;
        break;
    case 3:
        intraBC_ctrls->enabled             = pcs->sc_class1;
        intraBC_ctrls->ibc_shift           = 1;
        intraBC_ctrls->ibc_direction       = 1;
        intraBC_ctrls->hash_4x4_blocks     = !svt_aom_get_disallow_4x4(pcs->enc_mode, pcs->slice_type);
        intraBC_ctrls->max_block_size_hash = scs->super_block_size;
        break;
    case 4:
        intraBC_ctrls->enabled             = pcs->sc_class1;
        intraBC_ctrls->ibc_shift           = 1;
        intraBC_ctrls->ibc_direction       = 0;
        intraBC_ctrls->hash_4x4_blocks     = !svt_aom_get_disallow_4x4(pcs->enc_mode, pcs->slice_type);
        intraBC_ctrls->max_block_size_hash = block_size_wide[BLOCK_16X16];
        break;
    case 5:
        intraBC_ctrls->enabled             = pcs->sc_class1;
        intraBC_ctrls->ibc_shift           = 1;
        intraBC_ctrls->ibc_direction       = 0;
        intraBC_ctrls->hash_4x4_blocks     = !svt_aom_get_disallow_4x4(pcs->enc_mode, pcs->slice_type);
        intraBC_ctrls->max_block_size_hash = block_size_wide[BLOCK_8X8];
        break;
    case 6:
        intraBC_ctrls->enabled             = pcs->sc_class1;
        intraBC_ctrls->ibc_shift           = 1;
        intraBC_ctrls->ibc_direction       = 1;
        intraBC_ctrls->hash_4x4_blocks     = !svt_aom_get_disallow_4x4(pcs->enc_mode, pcs->slice_type);
        intraBC_ctrls->max_block_size_hash = block_size_wide[BLOCK_8X8];
        break;
    default: assert(0); break;
    }
}
/*
    set controls for Palette prediction
*/
static void set_palette_level(PictureParentControlSet *pcs, uint8_t palette_level) {
    PaletteCtrls *palette_ctrls = &pcs->palette_ctrls;

    switch (palette_level) {
    case 0: palette_ctrls->enabled = 0; break;
    case 1:
        palette_ctrls->enabled                       = 1;
        palette_ctrls->dominant_color_step           = 1;
        palette_ctrls->reduce_palette_cost_precision = 0;
        break;
    case 2:
        palette_ctrls->enabled                       = 1;
        palette_ctrls->dominant_color_step           = 2;
        palette_ctrls->reduce_palette_cost_precision = 0;
        break;
    case 3:
        palette_ctrls->enabled                       = 1;
        palette_ctrls->dominant_color_step           = 2;
        palette_ctrls->reduce_palette_cost_precision = 1;
        break;
    default: assert(0); break;
    }
}

/*
* return the max canidate count for MDS0
  Used by candidate injection and memory allocation
*/
uint16_t svt_aom_get_max_can_count(EncMode enc_mode) {
    //NOTE: this is a memory feature and not a speed feature. it should not be have any speed/quality impact.
    uint16_t mem_max_can_count;

    if (enc_mode <= ENC_M0)
        mem_max_can_count = 1225;
    else if (enc_mode <= ENC_M1)
        mem_max_can_count = 1000;
    else if (enc_mode <= ENC_M2)
        mem_max_can_count = 720;
    else if (enc_mode <= ENC_M3)
        mem_max_can_count = 576;
    else if (enc_mode <= ENC_M4)
        mem_max_can_count = 369;
    else if (enc_mode <= ENC_M5)
        mem_max_can_count = 295;
    else if (enc_mode <= ENC_M6)
        mem_max_can_count = 236;
    else if (enc_mode <= ENC_M11)
        mem_max_can_count = 190;
    else
        mem_max_can_count = 80;

    return mem_max_can_count;
}

/******************************************************
* Derive Multi-Processes Settings for OQ
Input   : encoder mode and tune
Output  : Multi-Processes signal(s)
******************************************************/
void svt_aom_sig_deriv_multi_processes(SequenceControlSet *scs, PictureParentControlSet *pcs,
                                       PictureDecisionContext *context_ptr) {
    FrameHeader            *frm_hdr             = &pcs->frm_hdr;
    EncMode                 enc_mode            = pcs->enc_mode;
    const uint8_t           is_islice           = pcs->slice_type == I_SLICE;
    const uint8_t           is_ref              = pcs->is_ref;
    const uint8_t           is_base             = pcs->temporal_layer_index == 0;
    const EbInputResolution input_resolution    = pcs->input_resolution;
    const uint32_t          bit_depth           = scs->static_config.encoder_bit_depth;
    const Bool              fast_decode         = scs->static_config.fast_decode;
    const uint32_t          hierarchical_levels = pcs->hierarchical_levels;
    const bool              rtc_tune  = (scs->static_config.pred_structure == SVT_AV1_PRED_LOW_DELAY_B) ? true : false;
    const uint8_t           sc_class1 = pcs->sc_class1;
    // Set GM ctrls assuming super-res is off for gm-pp need
    svt_aom_set_gm_controls(pcs, svt_aom_derive_gm_level(pcs, true));

    // If enabled here, the hme enable flags should also be enabled in ResourceCoordinationProcess
    // to ensure that resources are allocated for the downsampled pictures used in HME
    pcs->enable_hme_flag        = 1;
    pcs->enable_hme_level0_flag = 1;
    if (sc_class1) {
        pcs->enable_hme_level1_flag = 1;
        pcs->enable_hme_level2_flag = 1;
    } else if (enc_mode <= ENC_M7) {
        pcs->enable_hme_level1_flag = 1;
        pcs->enable_hme_level2_flag = 1;
    } else {
        pcs->enable_hme_level1_flag = 1;
        pcs->enable_hme_level2_flag = 0;
    }
    switch (pcs->tf_ctrls.hme_me_level) {
    case 0:
        pcs->tf_enable_hme_flag        = 1;
        pcs->tf_enable_hme_level0_flag = 1;
        pcs->tf_enable_hme_level1_flag = 1;
        pcs->tf_enable_hme_level2_flag = 1;
        break;

    case 1:
        pcs->tf_enable_hme_flag        = 1;
        pcs->tf_enable_hme_level0_flag = 1;
        pcs->tf_enable_hme_level1_flag = 1;
        pcs->tf_enable_hme_level2_flag = 0;
        break;

    case 2:
    case 3:
        pcs->tf_enable_hme_flag        = 1;
        pcs->tf_enable_hme_level0_flag = 1;
        pcs->tf_enable_hme_level1_flag = 0;
        pcs->tf_enable_hme_level2_flag = 0;
        break;

    default: assert(0); break;
    }
    // Set the Multi-Pass PD level
    pcs->multi_pass_pd_level = MULTI_PASS_PD_ON;
    // Loop filter Level                            Settings
    // 0                                            OFF
    // 1                                            CU-BASED
    // 2                                            LIGHT FRAME-BASED
    // 3                                            FULL FRAME-BASED

    //for now only I frames are allowed to use sc tools.
    //TODO: we can force all frames in GOP with the same detection status of leading I frame.
    uint8_t intrabc_level = 0;
    if (is_islice) {
        if (scs->intrabc_mode == DEFAULT) {
            if (rtc_tune) {
                intrabc_level = 0;
            } else {
                if (enc_mode <= ENC_M5)
                    intrabc_level = 1;
                else if (enc_mode <= ENC_M7)
                    intrabc_level = 4;
                else if (enc_mode <= ENC_M9)
                    intrabc_level = 5;
                else if (enc_mode <= ENC_M10)
                    intrabc_level = 6;
                else
                    intrabc_level = 0;
            }
        } else {
            intrabc_level = (uint8_t)scs->intrabc_mode;
        }
    } else {
        //this will enable sc tools for P frames. hence change Bitstream even if palette mode is OFF
        intrabc_level = 0;
    }
    set_intrabc_level(pcs, scs, intrabc_level);
    frm_hdr->allow_intrabc = pcs->intraBC_ctrls.enabled;
    // Set palette_level
    if (sc_class1) {
        if (scs->palette_level == DEFAULT) { //auto mode; if not set by cfg
            if (rtc_tune)
                pcs->palette_level = is_islice && sc_class1 ? 3 : is_islice ? 2 : 0;
            else if (enc_mode <= ENC_M11)
                pcs->palette_level = is_base ? 2 : 0;
            else if (enc_mode <= ENC_M12)
                pcs->palette_level = is_islice ? 2 : 0;
            else
                pcs->palette_level = 0;
        } else
            pcs->palette_level = scs->palette_level;
    } else
        pcs->palette_level = 0;

    set_palette_level(pcs, pcs->palette_level);
    uint8_t dlf_level = 0;
    if (pcs->scs->static_config.enable_dlf_flag && frm_hdr->allow_intrabc == 0) {
        dlf_level = get_dlf_level(
            enc_mode, is_ref, bit_depth > EB_EIGHT_BIT, fast_decode, input_resolution, rtc_tune, sc_class1);
    }

    frm_hdr->allow_screen_content_tools = sc_class1 && pcs->palette_level > 0 ? 1 : 0;

    svt_aom_set_dlf_controls(pcs, dlf_level);

    // Set CDEF controls
    if (scs->seq_header.cdef_level && frm_hdr->allow_intrabc == 0) {
        if (scs->static_config.cdef_level == DEFAULT) {
            if (rtc_tune) {
                if (enc_mode <= ENC_M11)
                    pcs->cdef_level = is_base ? 8 : is_ref ? 9 : 10;
                else {
                    if (input_resolution <= INPUT_SIZE_1080p_RANGE)
                        pcs->cdef_level = is_base ? 19 : 20;
                    else
                        pcs->cdef_level = is_islice ? 15 : is_ref ? 16 : 17;
                }
            } else if (enc_mode <= ENC_MR)
                pcs->cdef_level = 1;
            else if (enc_mode <= ENC_M2)
                pcs->cdef_level = 2;
            else if (enc_mode <= ENC_M5)
                pcs->cdef_level = 4;
            else if (enc_mode <= ENC_M9)
                pcs->cdef_level = is_base ? 8 : is_ref ? 9 : 10;
            else if (enc_mode <= ENC_M11)
                pcs->cdef_level = is_base ? 15 : is_ref ? 16 : 17;
            else if (enc_mode <= ENC_M13) {
                if (input_resolution <= INPUT_SIZE_1080p_RANGE)
                    pcs->cdef_level = is_base ? 15 : is_ref ? 16 : 17;
                else
                    pcs->cdef_level = is_islice ? 15 : is_ref ? 16 : 17;
            } else {
                if (hierarchical_levels <= 2) {
                    pcs->cdef_level = is_islice ? 15 : is_base ? 16 : 0;
                } else {
                    if (input_resolution <= INPUT_SIZE_1080p_RANGE)
                        pcs->cdef_level = is_base ? 15 : 0;
                    else
                        pcs->cdef_level = is_islice ? 15 : 0;
                }
            }
        } else
            pcs->cdef_level = (int8_t)(scs->static_config.cdef_level);
    } else
        pcs->cdef_level = 0;
    set_cdef_controls(pcs, pcs->cdef_level, fast_decode);

    uint8_t wn = 0, sg = 0;
    // If restoration filtering is enabled at the sequence level, derive the settings used for this frame
    if (scs->seq_header.enable_restoration) {
        wn = svt_aom_get_wn_filter_level(enc_mode, input_resolution, is_ref);
        sg = svt_aom_get_sg_filter_level(enc_mode, fast_decode, input_resolution, is_base);
    }

    Av1Common *cm = pcs->av1_cm;
    svt_aom_set_wn_filter_ctrls(cm, wn);
    svt_aom_set_sg_filter_ctrls(cm, sg);

    // Set whether restoration filtering is enabled for this frame
    pcs->enable_restoration = (wn > 0 || sg > 0);

    // Set frame end cdf update mode      Settings
    // 0                                     OFF
    // 1                                     ON
    if (scs->frame_end_cdf_update == DEFAULT)
        pcs->frame_end_cdf_update_mode = 1;
    else
        pcs->frame_end_cdf_update_mode = scs->frame_end_cdf_update;

    (void)context_ptr;

    // Tune TPL for better chroma.Only for 240P. 0 is OFF
#if TUNE_CHROMA_SSIM
    pcs->tune_tpl_for_chroma = 1;
#else
    pcs->tune_tpl_for_chroma = 0;
#endif
    if (scs->enable_hbd_mode_decision == DEFAULT)
        if (enc_mode <= ENC_M1)
            pcs->hbd_md = is_ref ? 1 : 2;
        else if (enc_mode <= ENC_M3)
            pcs->hbd_md = 2;
        else if (enc_mode <= ENC_M4)
            pcs->hbd_md = is_ref ? 2 : 0;
        else if (enc_mode <= ENC_M7)
            pcs->hbd_md = is_base ? 2 : 0;
        else
            pcs->hbd_md = is_islice ? 2 : 0;
    else
        pcs->hbd_md = scs->enable_hbd_mode_decision;

    pcs->max_can_count = svt_aom_get_max_can_count(enc_mode);
    if (enc_mode <= ENC_M8)
        pcs->use_best_me_unipred_cand_only = 0;
    else
        pcs->use_best_me_unipred_cand_only = 1;

    //TPL level should not be modified outside of this function
    svt_aom_set_tpl_extended_controls(pcs, scs->tpl_level);
    pcs->r0_based_qps_qpm = pcs->tpl_ctrls.enable &&
        (pcs->temporal_layer_index == 0 ||
         (scs->static_config.rate_control_mode == SVT_AV1_RC_MODE_CQP_OR_CRF &&
          ((pcs->hierarchical_levels == 5 && pcs->temporal_layer_index <= 2) ||
           (pcs->hierarchical_levels >= 4 && pcs->temporal_layer_index <= 1))));

    // If TPL results are needed for the current hierarchical layer, but are not available, shut r0-based QPS/QPM
    if (pcs->r0_based_qps_qpm && pcs->tpl_ctrls.reduced_tpl_group >= 0 &&
        pcs->temporal_layer_index > pcs->tpl_ctrls.reduced_tpl_group) {
        assert(pcs->temporal_layer_index != 0);
        pcs->r0_based_qps_qpm = 0;
    }
    pcs->adjust_under_shoot_gf = 0;
    if (scs->passes == 1 && scs->static_config.rate_control_mode == SVT_AV1_RC_MODE_VBR)
        pcs->adjust_under_shoot_gf = enc_mode <= ENC_M11 ? 1 : 2;
    const uint8_t ld_enhanced_base_frame_interval = 12;
    pcs->ld_enhanced_base_frame =
        (rtc_tune && is_base && !is_islice &&
         ((pcs->picture_number - scs->enc_ctx->last_idr_picture) % ld_enhanced_base_frame_interval == 0))
        ? 1
        : 0;
    if (enc_mode <= ENC_M12 || !rtc_tune)
        pcs->update_ref_count = 0;
    else
        pcs->update_ref_count = 1;
}
/*set the skip_frame flag for ipp*/
static void svt_aom_set_skip_frame_in_ipp(PictureParentControlSet *pcs, MeContext *me_ctx) {
    me_ctx->skip_frame             = 0;
    int8_t   max_heirachical_level = (int8_t)pcs->scs->max_heirachical_level;
    int8_t   enc_mode              = pcs->scs->static_config.enc_mode;
    int8_t   pass                  = (int8_t)pcs->scs->static_config.pass;
    uint8_t  skip_frame_first_pass = pcs->scs->ipp_pass_ctrls.skip_frame_first_pass;
    uint64_t picture_number        = pcs->picture_number;
    if (!skip_frame_first_pass) {
        if (enc_mode < ENC_M8)
            me_ctx->skip_frame = 0;
        else if (pass == ENC_SINGLE_PASS) {
            // SINGLE PASS
            if (max_heirachical_level <= 2)
                me_ctx->skip_frame = ((enc_mode >= ENC_M10) && (picture_number > 1 && picture_number % 2 > 0)) ? 1 : 0;
            else if ((enc_mode < ENC_M10 && (picture_number > 3 && picture_number % 4 > 0)) ||
                     (enc_mode >= ENC_M10 && (picture_number > 5 && picture_number % 6 > 0))) {
                me_ctx->skip_frame = 1;
            }
        }
    } else {
        if ((skip_frame_first_pass == 1) && (picture_number % 8 > 0))
            me_ctx->skip_frame = 1;
        else if ((skip_frame_first_pass == 2) && (picture_number > 7 && picture_number % 8 > 0))
            me_ctx->skip_frame = 1;
        else if (picture_number > 3 && picture_number % 4 > 0)
            me_ctx->skip_frame = 1;
        else
            me_ctx->skip_frame = 0;
    }
}

/******************************************************
* GM controls
******************************************************/
void svt_aom_set_gm_controls(PictureParentControlSet *pcs, uint8_t gm_level) {
    GmControls *gm_ctrls = &pcs->gm_ctrls;
    switch (gm_level) {
    case 0:
        gm_ctrls->enabled      = 0;
        gm_ctrls->use_ref_info = 0;
        gm_ctrls->pp_enabled   = 0;
        break;

    case 1:
        gm_ctrls->enabled                      = 1;
        gm_ctrls->identiy_exit                 = 0;
        gm_ctrls->rotzoom_model_only           = 0;
        gm_ctrls->bipred_only                  = 0;
        gm_ctrls->bypass_based_on_me           = 0;
        gm_ctrls->use_stationary_block         = 0;
        gm_ctrls->use_distance_based_active_th = 0;
        gm_ctrls->params_refinement_steps      = 5;
        gm_ctrls->downsample_level             = GM_FULL;

        gm_ctrls->corners       = 4;
        gm_ctrls->chess_rfn     = 0;
        gm_ctrls->match_sz      = 13;
        gm_ctrls->inj_psq_glb   = FALSE;
        gm_ctrls->use_ref_info  = 0;
        gm_ctrls->layer_offset  = 0;
        gm_ctrls->pp_enabled    = 0;
        gm_ctrls->ref_idx0_only = 0;
        break;
    case 2:
        gm_ctrls->enabled                      = 1;
        gm_ctrls->identiy_exit                 = 1;
        gm_ctrls->rotzoom_model_only           = 1;
        gm_ctrls->bipred_only                  = 0;
        gm_ctrls->bypass_based_on_me           = 0;
        gm_ctrls->use_stationary_block         = 0;
        gm_ctrls->use_distance_based_active_th = 0;
        gm_ctrls->params_refinement_steps      = 5;
        gm_ctrls->downsample_level             = GM_FULL;
        gm_ctrls->corners                      = 2;
        gm_ctrls->chess_rfn                    = 0;
        gm_ctrls->match_sz                     = 7;
        gm_ctrls->inj_psq_glb                  = FALSE;
        gm_ctrls->use_ref_info                 = 0;
        gm_ctrls->layer_offset                 = 0;
        gm_ctrls->pp_enabled                   = 1;
        gm_ctrls->ref_idx0_only                = 0;
        break;
    case 3:
        gm_ctrls->enabled                      = 1;
        gm_ctrls->identiy_exit                 = 1;
        gm_ctrls->rotzoom_model_only           = 1;
        gm_ctrls->bipred_only                  = 0;
        gm_ctrls->bypass_based_on_me           = 0;
        gm_ctrls->use_stationary_block         = 0;
        gm_ctrls->use_distance_based_active_th = 0;
        gm_ctrls->params_refinement_steps      = 5;
        gm_ctrls->downsample_level             = GM_FULL;
        gm_ctrls->corners                      = 2;
        gm_ctrls->chess_rfn                    = 0;
        gm_ctrls->match_sz                     = 7;
        gm_ctrls->inj_psq_glb                  = FALSE;
        gm_ctrls->use_ref_info                 = 0;
        gm_ctrls->layer_offset                 = 0;
        gm_ctrls->pp_enabled                   = 1;
        gm_ctrls->ref_idx0_only                = 0;
        break;
    case 4:
        gm_ctrls->enabled                      = 1;
        gm_ctrls->identiy_exit                 = 1;
        gm_ctrls->rotzoom_model_only           = 1;
        gm_ctrls->bipred_only                  = 0;
        gm_ctrls->bypass_based_on_me           = 1;
        gm_ctrls->use_stationary_block         = 0;
        gm_ctrls->use_distance_based_active_th = 0;
        gm_ctrls->params_refinement_steps      = 5;
        gm_ctrls->downsample_level             = GM_ADAPT_0;
        gm_ctrls->corners                      = 2;
        gm_ctrls->chess_rfn                    = 0;
        gm_ctrls->match_sz                     = 7;
        gm_ctrls->inj_psq_glb                  = TRUE;
        gm_ctrls->use_ref_info                 = 0;
        gm_ctrls->layer_offset                 = 0;
        gm_ctrls->pp_enabled                   = 1;
        gm_ctrls->ref_idx0_only                = 1;
        break;
    case 5:
        gm_ctrls->enabled                      = 1;
        gm_ctrls->identiy_exit                 = 1;
        gm_ctrls->rotzoom_model_only           = 1;
        gm_ctrls->bipred_only                  = 0;
        gm_ctrls->bypass_based_on_me           = 1;
        gm_ctrls->use_stationary_block         = 0;
        gm_ctrls->use_distance_based_active_th = 0;
        gm_ctrls->params_refinement_steps      = 5;
        gm_ctrls->downsample_level             = GM_ADAPT_1;
        gm_ctrls->corners                      = 2;
        gm_ctrls->chess_rfn                    = 1;
        gm_ctrls->match_sz                     = 7;
        gm_ctrls->inj_psq_glb                  = TRUE;
        gm_ctrls->use_ref_info                 = 0;
        gm_ctrls->layer_offset                 = 0;
        gm_ctrls->pp_enabled                   = 1;
        gm_ctrls->ref_idx0_only                = 1;
        break;
    case 6:
        gm_ctrls->enabled                      = 1;
        gm_ctrls->identiy_exit                 = 1;
        gm_ctrls->rotzoom_model_only           = 1;
        gm_ctrls->bipred_only                  = 0;
        gm_ctrls->bypass_based_on_me           = 1;
        gm_ctrls->use_stationary_block         = 0;
        gm_ctrls->use_distance_based_active_th = 0;
        gm_ctrls->params_refinement_steps      = 5;
        gm_ctrls->downsample_level             = GM_DOWN16;
        gm_ctrls->corners                      = 2;
        gm_ctrls->chess_rfn                    = 1;
        gm_ctrls->match_sz                     = 7;
        gm_ctrls->inj_psq_glb                  = TRUE;
        gm_ctrls->use_ref_info                 = 0;
        gm_ctrls->layer_offset                 = 3;
        gm_ctrls->pp_enabled                   = 0;
        gm_ctrls->ref_idx0_only                = 1;
        break;
    case 7:
        gm_ctrls->enabled                      = 1;
        gm_ctrls->identiy_exit                 = 1;
        gm_ctrls->rotzoom_model_only           = 1;
        gm_ctrls->bipred_only                  = 0;
        gm_ctrls->bypass_based_on_me           = 1;
        gm_ctrls->use_stationary_block         = 1;
        gm_ctrls->use_distance_based_active_th = 1;
        gm_ctrls->params_refinement_steps      = 5;
        gm_ctrls->downsample_level             = GM_DOWN16;
        gm_ctrls->corners                      = 2;
        gm_ctrls->chess_rfn                    = 1;
        gm_ctrls->match_sz                     = 7;
        gm_ctrls->inj_psq_glb                  = TRUE;
        gm_ctrls->use_ref_info                 = 0; //TODO: clean up ref info method
        gm_ctrls->layer_offset                 = 3;
        gm_ctrls->pp_enabled                   = 0;
        gm_ctrls->ref_idx0_only                = 1;
        break;
    case 8:
        gm_ctrls->enabled                      = 1;
        gm_ctrls->identiy_exit                 = 1;
        gm_ctrls->rotzoom_model_only           = 1;
        gm_ctrls->bipred_only                  = 1;
        gm_ctrls->bypass_based_on_me           = 1;
        gm_ctrls->use_stationary_block         = 1;
        gm_ctrls->use_distance_based_active_th = 1;
        gm_ctrls->params_refinement_steps      = 1;
        gm_ctrls->downsample_level             = GM_DOWN16;
        gm_ctrls->corners                      = 2;
        gm_ctrls->chess_rfn                    = 1;
        gm_ctrls->match_sz                     = 7;
        gm_ctrls->inj_psq_glb                  = TRUE;
        gm_ctrls->use_ref_info                 = 0;
        gm_ctrls->layer_offset                 = 3;
        gm_ctrls->pp_enabled                   = 0;
        gm_ctrls->ref_idx0_only                = 1;
        break;
    default: assert(0); break;
    }
    if (gm_level)
        assert((gm_ctrls->match_sz & 1) == 1);
}
/************************************************
 * Set ME/HME Params for the first pass encoding
 ************************************************/
static void *set_first_pass_me_hme_params_oq(MeContext *me_ctx, SequenceControlSet *scs,
                                             EbInputResolution input_resolution) {
    // HME/ME default settings
    me_ctx->num_hme_sa_w            = 2;
    me_ctx->num_hme_sa_h            = 2;
    me_ctx->reduce_hme_l0_sr_th_min = 0;
    me_ctx->reduce_hme_l0_sr_th_max = 0;

    // Set the minimum ME search area
    if (!scs->ipp_pass_ctrls.reduce_me_search) {
        me_ctx->me_sa.sa_min = (SearchArea){8, 8};
        me_ctx->me_sa.sa_max = (SearchArea){8, 8};
    } else {
        me_ctx->me_sa.sa_min = (SearchArea){8, 3};
        me_ctx->me_sa.sa_max = (SearchArea){8, 5};
    }
    if (input_resolution < INPUT_SIZE_1080p_RANGE) {
        me_ctx->hme_l0_sa.sa_min = (SearchArea){4, 4};
        me_ctx->hme_l0_sa.sa_max = (SearchArea){96, 96};
    } else {
        me_ctx->hme_l0_sa.sa_min = (SearchArea){8, 8};
        me_ctx->hme_l0_sa.sa_max = (SearchArea){96, 96};
    }
    me_ctx->hme_l1_sa = (SearchArea){8, 8};
    me_ctx->hme_l2_sa = (SearchArea){8, 8};
    // Scale up the MIN ME area if low frame rate
    bool low_frame_rate_flag = (scs->frame_rate >> 16);
    if (low_frame_rate_flag) {
        me_ctx->me_sa.sa_min.width  = (me_ctx->me_sa.sa_min.width * 3) >> 1;
        me_ctx->me_sa.sa_min.height = (me_ctx->me_sa.sa_min.height * 3) >> 1;
    }
    me_ctx->me_early_exit_th    = 0;
    me_ctx->me_safe_limit_zz_th = 0;
    return NULL;
};
/******************************************************
* Derive Pre-Analysis settings for first pass for pcs
Input   : encoder mode and tune
Output  : Pre-Analysis signal(s)
******************************************************/
void svt_aom_first_pass_sig_deriv_pre_analysis_pcs(PictureParentControlSet *pcs) {
    // Derive HME Flag
    pcs->enable_hme_flag        = 1;
    pcs->enable_hme_level0_flag = 1;
    pcs->enable_hme_level1_flag = 1;
    pcs->enable_hme_level2_flag = 1;

    //// Set here to allocate resources for the downsampled pictures used in HME (generated in PictureAnalysis)
    //// Will be later updated for SC/NSC in PictureDecisionProcess
    pcs->tf_enable_hme_flag        = 0;
    pcs->tf_enable_hme_level0_flag = 0;
    pcs->tf_enable_hme_level1_flag = 0;
    pcs->tf_enable_hme_level2_flag = 0;
}

/******************************************************
* Derive Pre-Analysis settings for first pass for scs
Input   : encoder mode and tune
Output  : Pre-Analysis signal(s)
******************************************************/
void svt_aom_first_pass_sig_deriv_pre_analysis_scs(SequenceControlSet *scs) {
    scs->seq_header.enable_intra_edge_filter   = 0;
    scs->seq_header.pic_based_rate_est         = 0;
    scs->seq_header.enable_restoration         = 0;
    scs->seq_header.cdef_level /*enable_cdef*/ = 0;
    scs->seq_header.enable_warped_motion       = 0;

    scs->seq_header.enable_superres                 = 0;
    scs->compound_mode                              = 0;
    scs->seq_header.order_hint_info.enable_jnt_comp = 0;
    scs->seq_header.enable_masked_compound          = 0;
    scs->seq_header.filter_intra_level              = 0;
    scs->seq_header.enable_interintra_compound      = 0;

    // Set the SCD Mode
    scs->scd_mode = scs->static_config.scene_change_detection == 0 ? SCD_MODE_0 : SCD_MODE_1;
}

/******************************************************
* Derive Multi-Processes Settings for first pass
Input   : encoder mode and tune
Output  : Multi-Processes signal(s)
******************************************************/

void svt_aom_first_pass_sig_deriv_multi_processes(SequenceControlSet *scs, PictureParentControlSet *pcs) {
    FrameHeader *frm_hdr = &pcs->frm_hdr;
    // If enabled here, the hme enable flags should also be enabled in ResourceCoordinationProcess
    // to ensure that resources are allocated for the downsampled pictures used in HME
    pcs->enable_hme_flag        = 1;
    pcs->enable_hme_level0_flag = 1;
    pcs->enable_hme_level1_flag = 1;
    pcs->enable_hme_level2_flag = 1;

    pcs->tf_enable_hme_flag        = 0;
    pcs->tf_enable_hme_level0_flag = 0;
    pcs->tf_enable_hme_level1_flag = 0;
    pcs->tf_enable_hme_level2_flag = 0;

    // Set the Multi-Pass PD level
    pcs->multi_pass_pd_level = MULTI_PASS_PD_OFF;

    // Set disallow_nsq

    pcs->max_number_of_pus_per_sb       = SQUARE_PU_COUNT;
    frm_hdr->allow_screen_content_tools = 0;
    frm_hdr->allow_intrabc              = 0;
    pcs->palette_level                  = 0;

    svt_aom_set_dlf_controls(pcs, 0);

    pcs->cdef_level = 0;

    Av1Common *cm = pcs->av1_cm;
    svt_aom_set_wn_filter_ctrls(cm, 0);
    svt_aom_set_sg_filter_ctrls(cm, 0);
    pcs->enable_restoration = 0;

    pcs->intra_pred_mode = 3;

    // Set frame end cdf update mode      Settings
    // 0                                  OFF
    // 1                                  ON
    if (scs->frame_end_cdf_update == DEFAULT)
        pcs->frame_end_cdf_update_mode = 0;
    else
        pcs->frame_end_cdf_update_mode = scs->frame_end_cdf_update;

    pcs->frm_hdr.use_ref_frame_mvs = 0;

    // GM off
    svt_aom_set_gm_controls(pcs, 0);
    pcs->ld_enhanced_base_frame = 0;
    pcs->update_ref_count       = 0;
}
/******************************************************
* Derive Mode Decision Config Settings for first pass
Input   : encoder mode and tune
Output  : EncDec Kernel signal(s)
******************************************************/
void svt_aom_first_pass_sig_deriv_mode_decision_config(PictureControlSet *pcs) {
    // CDF
    pcs->cdf_ctrl.enabled = pcs->cdf_ctrl.update_coef = 0;
    pcs->cdf_ctrl.update_mv = pcs->cdf_ctrl.update_se = 0;

    // Filter INTRA
    // pic_filter_intra_level specifies whether filter intra would be active
    // for a given picture.
    // pic_filter_intra_level | Settings
    // 0                      | OFF
    // 1                      | ON
    pcs->pic_filter_intra_level = 0;

    // High Precision
    FrameHeader *frm_hdr             = &pcs->ppcs->frm_hdr;
    frm_hdr->allow_high_precision_mv = 0;

    // Warped
    frm_hdr->allow_warped_motion       = 0;
    frm_hdr->is_motion_mode_switchable = frm_hdr->allow_warped_motion;

    // pic_obmc_level - pic_obmc_level is used to define md_pic_obmc_level.
    // The latter determines the OBMC settings in the function set_obmc_controls.
    // Please check the definitions of the flags/variables in the function
    // set_obmc_controls corresponding to the pic_obmc_level settings.
    //  pic_obmc_level  | Default Encoder Settings
    //         0        | OFF subject to possible constraints
    //       > 1        | Faster level subject to possible constraints
    pcs->ppcs->pic_obmc_level = 0;

    // Switchable Motion Mode
    frm_hdr->is_motion_mode_switchable = frm_hdr->is_motion_mode_switchable || pcs->ppcs->pic_obmc_level;

    // HBD Mode
    pcs->hbd_md                      = EB_8_BIT_MD; //first pass hard coded to 8bit
    pcs->ppcs->use_accurate_part_ctx = TRUE;
    pcs->ppcs->bypass_cost_table_gen = 0;
    pcs->approx_inter_rate           = 0;
}
/******************************************************
* Derive ME Settings for first pass
  Input   : encoder mode and tune
  Output  : ME Kernel signal(s)
******************************************************/
void svt_aom_first_pass_sig_deriv_me(SequenceControlSet *scs, PictureParentControlSet *pcs, MeContext *me_ctx) {
    // Set ME/HME search regions
    set_first_pass_me_hme_params_oq(me_ctx, scs, scs->input_resolution);

    // Set HME flags
    me_ctx->enable_hme_flag        = pcs->enable_hme_flag;
    me_ctx->enable_hme_level0_flag = pcs->enable_hme_level0_flag;
    me_ctx->enable_hme_level1_flag = pcs->enable_hme_level1_flag;
    me_ctx->enable_hme_level2_flag = !scs->ipp_pass_ctrls.reduce_me_search ? pcs->enable_hme_level2_flag : 0;

    // HME Search Method
    me_ctx->hme_search_method = SUB_SAD_SEARCH;

    // ME Search Method
    me_ctx->me_search_method = SUB_SAD_SEARCH;

    uint8_t gm_level = 0;
    svt_aom_set_gm_controls(pcs, gm_level);

    // Set pre-hme level (0-2)
    uint8_t prehme_level = 0;
    svt_aom_set_prehme_ctrls(me_ctx, prehme_level);

    // Set hme/me based reference pruning level (0-4)
    svt_aom_set_me_hme_ref_prune_ctrls(me_ctx, 0);

    // Set hme-based me sr adjustment level
    svt_aom_set_me_sr_adjustment_ctrls(me_ctx, 0);
    svt_aom_set_me_8x8_var_ctrls(me_ctx, 0);
    me_ctx->prune_me_candidates_th     = 0; // No impact on tf
    me_ctx->use_best_unipred_cand_only = 0; // No impact on tf
    svt_aom_set_prehme_ctrls(me_ctx, 0);
    svt_aom_set_skip_frame_in_ipp(pcs, me_ctx);
};
static void set_inter_comp_controls(ModeDecisionContext *ctx, PictureControlSet *pcs, uint8_t inter_comp_mode) {
    InterCompCtrls *inter_comp_ctrls = &ctx->inter_comp_ctrls;

    switch (inter_comp_mode) {
    case 0: //OFF
        inter_comp_ctrls->tot_comp_types = 1;
        break;
    case 1: //FULL
        inter_comp_ctrls->tot_comp_types = (svt_aom_get_wedge_params_bits(ctx->blk_geom->bsize) == 0) ? MD_COMP_WEDGE
                                                                                                      : MD_COMP_TYPES;
        inter_comp_ctrls->do_nearest_nearest  = 1;
        inter_comp_ctrls->do_near_near        = 1;
        inter_comp_ctrls->do_me               = 1;
        inter_comp_ctrls->do_pme              = 1;
        inter_comp_ctrls->do_nearest_near_new = 1;
        inter_comp_ctrls->do_3x3_bi           = 1;

        inter_comp_ctrls->skip_mvp_on_ref_info = 0;
        inter_comp_ctrls->use_rate             = 1;
        inter_comp_ctrls->pred0_to_pred1_mult  = 0;
        inter_comp_ctrls->no_dist              = 0;

        break;
    case 2:
        inter_comp_ctrls->tot_comp_types = (svt_aom_get_wedge_params_bits(ctx->blk_geom->bsize) == 0) ? MD_COMP_WEDGE
                                                                                                      : MD_COMP_TYPES;
        inter_comp_ctrls->do_nearest_nearest  = 1;
        inter_comp_ctrls->do_near_near        = 1;
        inter_comp_ctrls->do_me               = 1;
        inter_comp_ctrls->do_pme              = 1;
        inter_comp_ctrls->do_nearest_near_new = 1;
        inter_comp_ctrls->do_3x3_bi           = 0;

        inter_comp_ctrls->skip_mvp_on_ref_info = 0;
        inter_comp_ctrls->use_rate             = 0;
        inter_comp_ctrls->pred0_to_pred1_mult  = 0;
        inter_comp_ctrls->no_dist              = 0;

        break;

    case 3:
        inter_comp_ctrls->tot_comp_types = (svt_aom_get_wedge_params_bits(ctx->blk_geom->bsize) == 0) ? MD_COMP_WEDGE
                                                                                                      : MD_COMP_TYPES;
        inter_comp_ctrls->do_nearest_nearest  = 1;
        inter_comp_ctrls->do_near_near        = 1;
        inter_comp_ctrls->do_me               = 1;
        inter_comp_ctrls->do_pme              = 1;
        inter_comp_ctrls->do_nearest_near_new = 0;
        inter_comp_ctrls->do_3x3_bi           = 0;

        inter_comp_ctrls->skip_mvp_on_ref_info = 0;
        inter_comp_ctrls->use_rate             = 0;
        inter_comp_ctrls->pred0_to_pred1_mult  = 1;
        inter_comp_ctrls->no_dist              = 0;

        break;
    case 4:
        inter_comp_ctrls->tot_comp_types = (svt_aom_get_wedge_params_bits(ctx->blk_geom->bsize) == 0) ? MD_COMP_WEDGE
                                                                                                      : MD_COMP_TYPES;
        inter_comp_ctrls->do_nearest_nearest  = 1;
        inter_comp_ctrls->do_near_near        = 1;
        inter_comp_ctrls->do_me               = 0;
        inter_comp_ctrls->do_pme              = 0;
        inter_comp_ctrls->do_nearest_near_new = 0;
        inter_comp_ctrls->do_3x3_bi           = 0;

        inter_comp_ctrls->skip_mvp_on_ref_info = 1;
        inter_comp_ctrls->use_rate             = 0;
        inter_comp_ctrls->pred0_to_pred1_mult  = 4;
        inter_comp_ctrls->no_dist              = 0;

        break;
    case 5:

        if (ctx->blk_geom->shape != PART_N || pcs->ppcs->is_gm_on == 1)
            inter_comp_ctrls->tot_comp_types = 1;
        else
            inter_comp_ctrls->tot_comp_types = (svt_aom_get_wedge_params_bits(ctx->blk_geom->bsize) == 0)
                ? MD_COMP_WEDGE
                : MD_COMP_TYPES;

        inter_comp_ctrls->do_nearest_nearest  = 1;
        inter_comp_ctrls->do_near_near        = 1;
        inter_comp_ctrls->do_me               = 0;
        inter_comp_ctrls->do_pme              = 0;
        inter_comp_ctrls->do_nearest_near_new = 0;
        inter_comp_ctrls->do_3x3_bi           = 0;

        inter_comp_ctrls->skip_mvp_on_ref_info = 1;
        inter_comp_ctrls->use_rate             = 0;
        inter_comp_ctrls->pred0_to_pred1_mult  = 4;
        inter_comp_ctrls->no_dist              = 1;
        break;
    default: assert(0); break;
    }
}

/******************************************************
* Derive md Settings(feature signals) that could be
  changed  at the block level
******************************************************/

void svt_aom_sig_deriv_block(PictureControlSet *pcs, ModeDecisionContext *ctx) {
    ctx->spatial_sse_full_loop_level = ctx->spatial_sse_ctrls.spatial_sse_full_loop_level;
    // set compound_types_to_try
    set_inter_comp_controls(ctx, pcs, ctx->inter_compound_mode);
    svt_aom_set_nic_controls(ctx, ctx->pd_pass == PD_PASS_0 ? 17 : pcs->nic_level);
    if (pcs->ppcs->gm_ctrls.enabled && pcs->ppcs->gm_ctrls.inj_psq_glb && ctx->blk_geom->shape != PART_N)
        if (ctx->avail_blk_flag[ctx->blk_geom->sqi_mds])
            if (ctx->md_blk_arr_nsq[ctx->blk_geom->sqi_mds].pred_mode != GLOBAL_GLOBALMV &&
                ctx->md_blk_arr_nsq[ctx->blk_geom->sqi_mds].pred_mode != GLOBALMV) {
                ctx->params_status       = 1;
                ctx->global_mv_injection = 0;
            }
}

uint8_t svt_aom_get_enable_sg(EncMode enc_mode, uint8_t input_resolution, Bool fast_decode) {
    uint8_t sg = 0;
    for (int is_base = 0; is_base < 2; is_base++) {
        sg = svt_aom_get_sg_filter_level(enc_mode, fast_decode, input_resolution, is_base);
        if (sg)
            break;
    }
    return (sg > 0);
}

/*
* return true if restoration filtering is enabled; false otherwise
  Used by signal_derivation_pre_analysis_oq and memory allocation
*/
uint8_t svt_aom_get_enable_restoration(EncMode enc_mode, int8_t config_enable_restoration, uint8_t input_resolution,
                                       Bool fast_decode) {
    if (config_enable_restoration != DEFAULT)
        return config_enable_restoration;

    uint8_t wn = 0;
    for (int is_ref = 0; is_ref < 2; is_ref++) {
        wn = svt_aom_get_wn_filter_level(enc_mode, input_resolution, is_ref);
        if (wn)
            break;
    }

    uint8_t sg = svt_aom_get_enable_sg(enc_mode, input_resolution, fast_decode);

    return (sg > 0 || wn > 0);
}

/******************************************************
* Derive Pre-Analysis settings for OQ for pcs
Input   : encoder mode and tune
Output  : Pre-Analysis signal(s)
******************************************************/
void svt_aom_sig_deriv_pre_analysis_pcs(PictureParentControlSet *pcs) {
    // Derive HME Flag
    // Set here to allocate resources for the downsampled pictures used in HME (generated in PictureAnalysis)
    // Will be later updated for SC/NSC in PictureDecisionProcess
    pcs->enable_hme_flag        = 1;
    pcs->enable_hme_level0_flag = 1;
    pcs->enable_hme_level1_flag = 1;
    pcs->enable_hme_level2_flag = 1;
    // Set here to allocate resources for the downsampled pictures used in HME (generated in PictureAnalysis)
    // Will be later updated for SC/NSC in PictureDecisionProcess
    pcs->tf_enable_hme_flag        = 1;
    pcs->tf_enable_hme_level0_flag = 1;
    pcs->tf_enable_hme_level1_flag = 1;
    pcs->tf_enable_hme_level2_flag = 1;

    //if (scs->static_config.enable_tpl_la)
    //svt_aom_assert_err(pcs->is_720p_or_larger == (pcs->tpl_ctrls.synth_blk_size == 16), "TPL Synth Size Error");
}

/******************************************************
* Derive Pre-Analysis settings for OQ for scs
Input   : encoder mode and tune
Output  : Pre-Analysis signal(s)
******************************************************/
void svt_aom_sig_deriv_pre_analysis_scs(SequenceControlSet *scs) {
    const int8_t enc_mode = scs->static_config.enc_mode;
    // Set the SCD Mode
    scs->scd_mode = scs->static_config.scene_change_detection == 0 ? SCD_MODE_0 : SCD_MODE_1;

    // initialize sequence level enable_superres
    scs->seq_header.enable_superres = scs->static_config.superres_mode > SUPERRES_NONE ? 1 : 0;
    uint8_t ii_allowed              = 0;
    for (uint8_t is_base = 0; is_base < 2; is_base++) {
        for (uint8_t transition_present = 0; transition_present < 2; transition_present++) {
            if (ii_allowed)
                break;
            ii_allowed |= svt_aom_get_inter_intra_level(enc_mode, is_base, transition_present);
        }
    }
    scs->seq_header.enable_interintra_compound = ii_allowed ? 1 : 0;

    if (get_filter_intra_level(enc_mode))
        scs->seq_header.filter_intra_level = 1;
    else
        scs->seq_header.filter_intra_level = 0;
    if (get_inter_compound_level(enc_mode)) {
        scs->seq_header.order_hint_info.enable_jnt_comp = 1; //DISTANCE
        scs->seq_header.enable_masked_compound          = 1; //DIFF+WEDGE
    } else {
        scs->seq_header.order_hint_info.enable_jnt_comp = 0;
        scs->seq_header.enable_masked_compound          = 0;
    }
    if (scs->enable_intra_edge_filter == DEFAULT)
        scs->seq_header.enable_intra_edge_filter = 1;
    else
        scs->seq_header.enable_intra_edge_filter = (uint8_t)scs->enable_intra_edge_filter;

    if (scs->pic_based_rate_est == DEFAULT)
        scs->seq_header.pic_based_rate_est = 0;
    else
        scs->seq_header.pic_based_rate_est = (uint8_t)scs->pic_based_rate_est;

    if (scs->static_config.enable_restoration_filtering == DEFAULT) {
        scs->seq_header.enable_restoration = svt_aom_get_enable_restoration(
            scs->static_config.enc_mode,
            scs->static_config.enable_restoration_filtering,
            scs->input_resolution,
            scs->static_config.fast_decode);
    } else
        scs->seq_header.enable_restoration = (uint8_t)scs->static_config.enable_restoration_filtering;

    if (scs->static_config.cdef_level == DEFAULT)
        scs->seq_header.cdef_level = 1;
    else
        scs->seq_header.cdef_level = (uint8_t)(scs->static_config.cdef_level > 0);

    if (scs->enable_warped_motion == DEFAULT) {
        scs->seq_header.enable_warped_motion = 1;
    } else
        scs->seq_header.enable_warped_motion = (uint8_t)scs->enable_warped_motion;
}

static uint32_t hadamard_path(ModeDecisionContext *ctx, uint8_t *input, uint32_t input_origin_index,
                              uint32_t input_stride, uint8_t *pred, uint32_t blk_origin_index, uint32_t pred_stride) {
    uint32_t input_idx, pred_idx, res_idx;

    uint32_t satd_cost = 0;

    const TxSize tx_size = AOMMIN(TX_32X32, max_txsize_lookup[ctx->blk_geom->bsize]);

    const int stepr = tx_size_high_unit[tx_size];
    const int stepc = tx_size_wide_unit[tx_size];
    const int txbw  = tx_size_wide[tx_size];
    const int txbh  = tx_size_high[tx_size];

    const int max_blocks_wide = block_size_wide[ctx->blk_geom->bsize] >> MI_SIZE_LOG2;
    const int max_blocks_high = block_size_wide[ctx->blk_geom->bsize] >> MI_SIZE_LOG2;
    int       row, col;

    for (row = 0; row < max_blocks_high; row += stepr) {
        for (col = 0; col < max_blocks_wide; col += stepc) {
            input_idx = input_origin_index + (((row * input_stride) + col) << 2);
            pred_idx  = blk_origin_index + (((row * pred_stride) + col) << 2);
            res_idx   = 0;

            svt_aom_residual_kernel(input,
                                    input_idx,
                                    input_stride,
                                    pred,
                                    pred_idx,
                                    pred_stride,
                                    (int16_t *)ctx->temp_residual->buffer_y,
                                    res_idx,
                                    ctx->temp_residual->stride_y,
                                    0, // input and pred 8-bit
                                    txbw,
                                    txbh);

            switch (tx_size) {
            case TX_4X4:
                svt_aom_hadamard_4x4((int16_t *)ctx->temp_residual->buffer_y,
                                     ctx->temp_residual->stride_y,
                                     &(((int32_t *)ctx->tx_coeffs->buffer_y)[0]));
                break;

            case TX_8X8:
                svt_aom_hadamard_8x8((int16_t *)ctx->temp_residual->buffer_y,
                                     ctx->temp_residual->stride_y,
                                     &(((int32_t *)ctx->tx_coeffs->buffer_y)[0]));
                break;

            case TX_16X16:
                svt_aom_hadamard_16x16((int16_t *)ctx->temp_residual->buffer_y,
                                       ctx->temp_residual->stride_y,
                                       &(((int32_t *)ctx->tx_coeffs->buffer_y)[0]));
                break;

            case TX_32X32:
                svt_aom_hadamard_32x32((int16_t *)ctx->temp_residual->buffer_y,
                                       ctx->temp_residual->stride_y,
                                       &(((int32_t *)ctx->tx_coeffs->buffer_y)[0]));
                break;

            default: assert(0);
            }
            satd_cost += svt_aom_satd(&(((int32_t *)ctx->tx_coeffs->buffer_y)[0]), tx_size_2d[tx_size]);
        }
    }
    return (satd_cost);
}

static EbErrorType svt_aom_check_high_freq(PictureControlSet *pcs, SuperBlock *sb_ptr, ModeDecisionContext *ctx) {
    EbErrorType         return_error   = EB_ErrorNone;
    SequenceControlSet *scs            = pcs->scs;
    uint8_t             ref_frame_type = NONE_FRAME;

    ctx->sb_ptr = sb_ptr;

    if (pcs->ppcs->me_32x32_distortion[ctx->sb_index] == 0 ||
        pcs->ppcs->me_8x8_cost_variance[ctx->sb_index] < ctx->detect_high_freq_ctrls.me_8x8_sad_var_th)
        return return_error;

    EbPictureBufferDesc *input_pic = pcs->ppcs->enhanced_pic;

    MotionEstimationData *pa_me_data                     = pcs->ppcs->pa_me_data;
    uint32_t              blk32_idx_tab[GEOM_TOT - 1][4] = {{1, 22, 43, 64},
                                                            {5, 30, 55, 80},
                                                            {5, 46, 87, 128},
                                                            {5, 110, 215, 320},
                                                            {5, 174, 343, 512},
                                                            {25, 294, 563, 832}};
    uint32_t              sum_b32_satd                   = 0;
    uint8_t               is_high_satd                   = 0; // the b64 is detected if at least one b32 is detected
    uint32_t              b32_satd[4];

    for (uint32_t blk_idx = 0; blk_idx < 4; blk_idx++) {
        b32_satd[blk_idx] = (uint32_t)~0;

        uint32_t blk_idx_mds = blk32_idx_tab[scs->svt_aom_geom_idx][blk_idx];

        // block position should be calculated from the values in MD context,
        // because sb params are different since frames might be downscaled
        // if super-res or resize is enabled
        const BlockGeom *blk_geom = ctx->blk_geom = get_blk_geom_mds(blk_idx_mds);
        ctx->blk_org_x                            = (uint16_t)(ctx->sb_origin_x + blk_geom->org_x);
        ctx->blk_org_y                            = (uint16_t)(ctx->sb_origin_y + blk_geom->org_y);
        const uint32_t input_origin_index         = (ctx->blk_org_y + input_pic->org_y) * input_pic->stride_y +
            (ctx->blk_org_x + input_pic->org_x);

        ctx->me_sb_addr = ctx->sb_ptr->index;

        ctx->me_block_offset = (ctx->blk_geom->svt_aom_geom_idx == GEOM_0) ? me_idx_85[ctx->blk_geom->blkidx_mds]
            : (ctx->blk_geom->svt_aom_geom_idx == GEOM_1)                  ? me_idx_geom1[ctx->blk_geom->blkidx_mds]
            : (ctx->blk_geom->svt_aom_geom_idx == GEOM_2)                  ? me_idx_geom2[ctx->blk_geom->blkidx_mds]
            : (ctx->blk_geom->svt_aom_geom_idx == GEOM_3)                  ? me_idx_geom3[ctx->blk_geom->blkidx_mds]
            : (ctx->blk_geom->svt_aom_geom_idx == GEOM_4)                  ? me_idx_geom4[ctx->blk_geom->blkidx_mds]
                                                                           : me_idx[ctx->blk_geom->blkidx_mds];

        ctx->me_cand_offset = ctx->me_block_offset * pa_me_data->max_cand;

        // ME offset(s)
        FrameHeader       *frm_hdr             = &pcs->ppcs->frm_hdr;
        Bool               is_compound_enabled = (frm_hdr->reference_mode == SINGLE_REFERENCE) ? 0 : 1;
        const MeSbResults *me_results          = pa_me_data->me_results[ctx->me_sb_addr];
        uint8_t            total_me_cnt        = me_results->total_me_candidate_index[ctx->me_block_offset];

        const MeCandidate *me_block_results = &me_results->me_candidate_array[ctx->me_cand_offset];

        const uint8_t max_refs = pa_me_data->max_refs;
        const uint8_t max_l0   = pa_me_data->max_l0;

        for (uint8_t me_candidate_index = 0; me_candidate_index < total_me_cnt; ++me_candidate_index) {
            const MeCandidate *me_block_results_ptr = &me_block_results[me_candidate_index];
            const uint8_t      inter_direction      = me_block_results_ptr->direction;

            if (inter_direction == 2)
                break;

            const uint8_t list0_ref_index = me_block_results_ptr->ref_idx_l0;
            const uint8_t list1_ref_index = me_block_results_ptr->ref_idx_l1;
            Mv            mv[MAX_NUM_OF_REF_PIC_LIST];
            /**************
                NEWMV L0
            ************* */
            mv[REF_LIST_0] = (Mv){{0, 0}};
            if (inter_direction == 0) {
                const int16_t mv_x = (me_results->me_mv_array[ctx->me_block_offset * max_refs + list0_ref_index].x_mv)
                    << 1;
                const int16_t mv_y = (me_results->me_mv_array[ctx->me_block_offset * max_refs + list0_ref_index].y_mv)
                    << 1;
                mv[REF_LIST_0] = (Mv){{mv_x, mv_y}};
                ref_frame_type = svt_get_ref_frame_type(REF_LIST_0, list0_ref_index);
            }

            /**************
                NEWMV L1
            ************* */
            mv[REF_LIST_1] = (Mv){{0, 0}};
            if (is_compound_enabled) {
                if (inter_direction == 1) {
                    const int16_t mv_x =
                        (me_results->me_mv_array[ctx->me_block_offset * max_refs + max_l0 + list0_ref_index].x_mv) << 1;
                    const int16_t mv_y =
                        (me_results->me_mv_array[ctx->me_block_offset * max_refs + max_l0 + list0_ref_index].y_mv) << 1;
                    mv[REF_LIST_1] = (Mv){{mv_x, mv_y}};
                    ref_frame_type = svt_get_ref_frame_type(REF_LIST_1, list1_ref_index);
                }
            }

            MvReferenceFrame rf[2];
            av1_set_ref_frame(rf, ref_frame_type);

            assert(rf[1] == NONE_FRAME);

            EbPictureBufferDesc *ref_pic;
            const int8_t         ref_idx_first  = get_ref_frame_idx(rf[0]);
            const int8_t         list_idx_first = get_list_idx(rf[0]);
            int32_t              ref_origin_index;
            int16_t              mv_x, mv_y;

            EbReferenceObject *ref_obj =
                (EbReferenceObject *)pcs->ref_pic_ptr_array[list_idx_first][ref_idx_first]->object_ptr;

            if (list_idx_first == 0) {
                ref_pic = svt_aom_get_ref_pic_buffer(pcs, 0, 0, ref_idx_first);
                mv_x    = mv[REF_LIST_0].x >> 3;
                mv_y    = mv[REF_LIST_0].y >> 3;
            } else {
                ref_pic = svt_aom_get_ref_pic_buffer(pcs, 0, 1, ref_idx_first);
                mv_x    = mv[REF_LIST_1].x >> 3;
                mv_y    = mv[REF_LIST_1].y >> 3;
            }
            // -------
            // Use scaled references if resolution of the reference is different from that of the input
            // -------
            svt_aom_use_scaled_rec_refs_if_needed(pcs, input_pic, ref_obj, &ref_pic, 0);
            ref_origin_index = ref_pic->org_x + (ctx->blk_org_x + mv_x) +
                (ctx->blk_org_y + mv_y + ref_pic->org_y) * ref_pic->stride_y;

            uint32_t satd = hadamard_path(ctx,
                                          input_pic->buffer_y,
                                          input_origin_index,
                                          input_pic->stride_y,
                                          ref_pic->buffer_y,
                                          ref_origin_index,
                                          ref_pic->stride_y);

            b32_satd[blk_idx] = MIN(b32_satd[blk_idx], satd);
        }

        if (b32_satd[blk_idx] >= ctx->detect_high_freq_ctrls.high_satd_th) {
            is_high_satd = 1;
        }
        sum_b32_satd += b32_satd[blk_idx];

        if (is_high_satd && sum_b32_satd > pcs->ppcs->me_32x32_distortion[ctx->sb_index]) {
            int dev = ((sum_b32_satd - pcs->ppcs->me_32x32_distortion[ctx->sb_index]) * 100) /
                pcs->ppcs->me_32x32_distortion[ctx->sb_index];
            if (dev >= ctx->detect_high_freq_ctrls.satd_to_sad_dev_th) {
                ctx->high_freq_present = 1;
                return return_error;
            }
        }
    }

    return return_error;
}

/*
* check if the reference picture is in same frame size
* TRUE -- in same frame size
* FALSE -- reference picture not exist or in difference frame size
*/
Bool svt_aom_is_ref_same_size(PictureControlSet *pcs, uint8_t list_idx, uint8_t ref_idx) {
    // skip the checking if reference scaling and super-res are disabled
    if (pcs->ppcs->is_not_scaled)
        return TRUE;
    if (pcs->slice_type != P_SLICE && pcs->slice_type != B_SLICE)
        return FALSE;
    if (pcs->slice_type != B_SLICE && list_idx == REF_LIST_1)
        return FALSE;
    if (pcs->ref_pic_ptr_array[list_idx][ref_idx] == NULL)
        return FALSE;

    EbReferenceObject *ref_obj = (EbReferenceObject *)pcs->ref_pic_ptr_array[list_idx][ref_idx]->object_ptr;
    if (ref_obj == NULL || ref_obj->reference_picture == NULL)
        return FALSE;

    return ref_obj->reference_picture->width == pcs->ppcs->frame_width &&
        ref_obj->reference_picture->height == pcs->ppcs->frame_height;
}
void svt_aom_set_obmc_controls(ModeDecisionContext *ctx, uint8_t obmc_mode) {
    ObmcControls *obmc_ctrls = &ctx->obmc_ctrls;
    switch (obmc_mode) {
    case 0: obmc_ctrls->enabled = 0; break;
    case 1:
        obmc_ctrls->enabled                      = 1;
        obmc_ctrls->max_blk_size_to_refine_16x16 = 0;
        obmc_ctrls->max_blk_size_16x16           = 0;
        obmc_ctrls->refine_level                 = 0;
        obmc_ctrls->trans_face_off               = 0;
        break;
    case 2:
        obmc_ctrls->enabled                      = 1;
        obmc_ctrls->max_blk_size_to_refine_16x16 = 1;
        obmc_ctrls->max_blk_size_16x16           = 0;
        obmc_ctrls->refine_level                 = 1;
        obmc_ctrls->trans_face_off               = 0;
        break;
    case 3:
        obmc_ctrls->enabled                      = 1;
        obmc_ctrls->max_blk_size_to_refine_16x16 = 1;
        obmc_ctrls->max_blk_size_16x16           = 0;
        obmc_ctrls->refine_level                 = 1;
        obmc_ctrls->trans_face_off               = 1;
        obmc_ctrls->trans_face_off_th            = 0;
        break;
    case 4:
        obmc_ctrls->enabled                      = 1;
        obmc_ctrls->max_blk_size_to_refine_16x16 = 1;
        obmc_ctrls->max_blk_size_16x16           = 1;
        obmc_ctrls->refine_level                 = 4;
        obmc_ctrls->trans_face_off               = 1;
        obmc_ctrls->trans_face_off_th            = 50;
        break;
    default: assert(0); break;
    }
}
void set_depth_removal_level_controls_rtc(PictureControlSet *pcs, ModeDecisionContext *ctx) {
    DepthRemovalCtrls *depth_removal_ctrls = &ctx->depth_removal_ctrls;
    B64Geom           *b64_geom            = &pcs->ppcs->b64_geom[ctx->sb_index];
    uint8_t            drl                 = pcs->pic_depth_removal_level_rtc;
    if (pcs->slice_type == I_SLICE) {
        switch (drl) {
        case 0: depth_removal_ctrls->enabled = 0; break;
        case 1:
            depth_removal_ctrls->enabled              = 1;
            depth_removal_ctrls->disallow_below_64x64 = 0;
            depth_removal_ctrls->disallow_below_32x32 = 0;
            depth_removal_ctrls->disallow_below_16x16 = 1;
            break;
        }
    } else {
        switch (drl) {
        case 0: depth_removal_ctrls->enabled = 0; break;

        case 1:
            depth_removal_ctrls->enabled              = 1;
            depth_removal_ctrls->disallow_below_64x64 = 0;
            depth_removal_ctrls->disallow_below_32x32 = 0;
            depth_removal_ctrls->disallow_below_16x16 = 0;
            {
                ctx->depth_removal_ctrls.enabled = 1;
                uint32_t var                     = pcs->ppcs->me_8x8_cost_variance[ctx->sb_index];
                uint32_t sad                     = pcs->ppcs->me_64x64_distortion[ctx->sb_index];

                //uint32_t var_th1 = 6 * ctx->qp_index * (pcs->temporal_layer_index + 1) * (1+pcs->parent_pcs_ptr->input_resolution);
                uint32_t var_th1 = ctx->qp_index * (pcs->temporal_layer_index + 1) * (7 + pcs->ppcs->input_resolution);
                uint32_t sad_th2 = var_th1 * 2;

                if (var < var_th1 && sad < sad_th2) {
                    ctx->depth_removal_ctrls.disallow_below_64x64 = 1;
                    ctx->depth_removal_ctrls.disallow_below_32x32 = 1;
                    ctx->depth_removal_ctrls.disallow_below_16x16 = 1;
                    ctx->skip_pd0                                 = 1;

                } else if (var < (var_th1 >> 3) && sad < (sad_th2 >> 3)) {
                    ctx->depth_removal_ctrls.disallow_below_64x64 = 0;
                    ctx->depth_removal_ctrls.disallow_below_32x32 = 1;
                    ctx->depth_removal_ctrls.disallow_below_16x16 = 1;
                    ctx->skip_pd0                                 = 0;
                } else {
                    ctx->depth_removal_ctrls.disallow_below_64x64 = 0;
                    ctx->depth_removal_ctrls.disallow_below_32x32 = 0;
                    ctx->depth_removal_ctrls.disallow_below_16x16 = 1;
                    ctx->skip_pd0                                 = 0;
                }
            }
            break;
        }
    }

    depth_removal_ctrls->disallow_below_16x16 = (b64_geom->width % 16 == 0 && b64_geom->height % 16 == 0)
        ? depth_removal_ctrls->disallow_below_16x16
        : 0;
    depth_removal_ctrls->disallow_below_32x32 = (b64_geom->width % 32 == 0 && b64_geom->height % 32 == 0)
        ? depth_removal_ctrls->disallow_below_32x32
        : 0;
    depth_removal_ctrls->disallow_below_64x64 = (b64_geom->width % 64 == 0 && b64_geom->height % 64 == 0)
        ? depth_removal_ctrls->disallow_below_64x64
        : 0;
}
static void set_depth_removal_level_controls(PictureControlSet *pcs, ModeDecisionContext *ctx,
                                             uint8_t depth_removal_level) {
    DepthRemovalCtrls *depth_removal_ctrls = &ctx->depth_removal_ctrls;

    if (pcs->slice_type == I_SLICE) {
        // Use b64_geom here because feature is only enabled when 64x64 SB size is used
        B64Geom *b64_geom = &pcs->ppcs->b64_geom[ctx->sb_index];

        uint16_t disallow_below_16x16_variance_th = 0;
        uint16_t disallow_below_32x32_variance_th = 0;
        uint16_t disallow_below_64x64_variance_th = 0;

        switch (depth_removal_level) {
        case 0: depth_removal_ctrls->enabled = 0; break;

        case 1:
            depth_removal_ctrls->enabled     = 1;
            disallow_below_16x16_variance_th = 150;
            disallow_below_32x32_variance_th = 50;
            disallow_below_64x64_variance_th = 25;
            break;
        }

        if (depth_removal_ctrls->enabled) {
            // If variance is available, use in depth removal decision
            if (pcs->ppcs->variance) {
                depth_removal_ctrls->disallow_below_16x16 = (b64_geom->width % 16 == 0 && b64_geom->height % 16 == 0)
                    ? (depth_removal_ctrls->disallow_below_16x16 ||
                       pcs->ppcs->variance[ctx->sb_index][ME_TIER_ZERO_PU_64x64] < disallow_below_16x16_variance_th)
                    : 0;

                depth_removal_ctrls->disallow_below_32x32 = (b64_geom->width % 32 == 0 && b64_geom->height % 32 == 0)
                    ? (depth_removal_ctrls->disallow_below_32x32 ||
                       pcs->ppcs->variance[ctx->sb_index][ME_TIER_ZERO_PU_64x64] < disallow_below_32x32_variance_th)
                    : 0;

                depth_removal_ctrls->disallow_below_64x64 = (b64_geom->width % 64 == 0 && b64_geom->height % 64 == 0)
                    ? (depth_removal_ctrls->disallow_below_64x64 ||
                       pcs->ppcs->variance[ctx->sb_index][ME_TIER_ZERO_PU_64x64] < disallow_below_64x64_variance_th)
                    : 0;
            } else {
                depth_removal_ctrls->disallow_below_16x16 = (b64_geom->width % 16 == 0 && b64_geom->height % 16 == 0)
                    ? depth_removal_ctrls->disallow_below_16x16
                    : 0;

                depth_removal_ctrls->disallow_below_32x32 = (b64_geom->width % 32 == 0 && b64_geom->height % 32 == 0)
                    ? depth_removal_ctrls->disallow_below_32x32
                    : 0;

                depth_removal_ctrls->disallow_below_64x64 = (b64_geom->width % 64 == 0 && b64_geom->height % 64 == 0)
                    ? depth_removal_ctrls->disallow_below_64x64
                    : 0;
            }
        }
    } else {
        uint32_t me_8x8_cost_variance = pcs->ppcs->me_8x8_cost_variance[ctx->sb_index];

        // Use b64_geom here because feature is only enabled when 64x64 SB size is used
        B64Geom *b64_geom = &pcs->ppcs->b64_geom[ctx->sb_index];

        // me_distortion => EB_8_BIT_MD
        uint32_t fast_lambda = ctx->fast_lambda_md[EB_8_BIT_MD];

        uint32_t sb_size = 64 * 64;

        uint64_t cost_th_rate = 1 << 13;

        uint64_t disallow_below_16x16_cost_th_multiplier = 0;
        uint64_t disallow_below_32x32_cost_th_multiplier = 0;
        uint64_t disallow_below_64x64_cost_th_multiplier = 0;

        int64_t dev_16x16_to_8x8_th   = MAX_SIGNED_VALUE;
        int64_t dev_32x32_to_16x16_th = MAX_SIGNED_VALUE;
        int64_t dev_32x32_to_8x8_th   = MAX_SIGNED_VALUE;

        int8_t qp_scale_factor = 0;

        // Modulate depth_removal level for Layer0 frames based on the qp_offset band
        if (pcs->ppcs->frm_hdr.delta_q_params.delta_q_present) {
            int diff = ctx->sb_ptr->qindex - quantizer_to_qindex[pcs->ppcs->picture_qp];
            if (diff <= -12)
                depth_removal_level = MAX(0, (int)depth_removal_level - 4);
            else if (diff <= -6)
                depth_removal_level = MAX(0, (int)depth_removal_level - 3);
            else if (diff <= -3)
                depth_removal_level = MAX(0, (int)depth_removal_level - 2);
            else if (diff < 0)
                depth_removal_level = MAX(0, (int)depth_removal_level - 1);
        }

        switch (depth_removal_level) {
        case 0: depth_removal_ctrls->enabled = 0; break;

        case 1:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 0;
            disallow_below_32x32_cost_th_multiplier = 0;
            disallow_below_64x64_cost_th_multiplier = 0;
            dev_16x16_to_8x8_th                     = 0;
            dev_32x32_to_16x16_th                   = 0;
            qp_scale_factor                         = 1;

            break;

        case 2:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 4;
            disallow_below_32x32_cost_th_multiplier = 0;
            disallow_below_64x64_cost_th_multiplier = 0;
            dev_16x16_to_8x8_th                     = 20;
            dev_32x32_to_16x16_th                   = 0;
            qp_scale_factor                         = 1;

            break;
        case 3:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 4;
            disallow_below_32x32_cost_th_multiplier = 0;
            disallow_below_64x64_cost_th_multiplier = 0;
            dev_16x16_to_8x8_th                     = 50;
            dev_32x32_to_16x16_th                   = 0;
            qp_scale_factor                         = 1;
            break;
        case 4:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 4;
            disallow_below_32x32_cost_th_multiplier = 0;
            disallow_below_64x64_cost_th_multiplier = 0;
            dev_16x16_to_8x8_th                     = 50;
            dev_32x32_to_16x16_th                   = 0;
            qp_scale_factor                         = 2;
            break;
        case 5:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 4;
            disallow_below_32x32_cost_th_multiplier = 0;
            disallow_below_64x64_cost_th_multiplier = 0;
            dev_16x16_to_8x8_th                     = 50;
            dev_32x32_to_16x16_th                   = 0;
            qp_scale_factor                         = 2;
            break;
        case 6:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 4;
            disallow_below_32x32_cost_th_multiplier = 2;
            disallow_below_64x64_cost_th_multiplier = 0;
            dev_16x16_to_8x8_th                     = 50;
            dev_32x32_to_16x16_th                   = 0;
            qp_scale_factor                         = 2;
            break;
        case 7:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 4;
            disallow_below_32x32_cost_th_multiplier = 2;
            disallow_below_64x64_cost_th_multiplier = 2;
            dev_16x16_to_8x8_th                     = 50;
            dev_32x32_to_16x16_th                   = 0;
            qp_scale_factor                         = 2;
            break;
        case 8:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 4;
            disallow_below_32x32_cost_th_multiplier = 2;
            disallow_below_64x64_cost_th_multiplier = 2;
            dev_16x16_to_8x8_th                     = 100;
            dev_32x32_to_16x16_th                   = 50;
            qp_scale_factor                         = 3;
            break;
        case 9:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 8;
            disallow_below_32x32_cost_th_multiplier = 2;
            disallow_below_64x64_cost_th_multiplier = 2;
            dev_16x16_to_8x8_th                     = 100;
            dev_32x32_to_16x16_th                   = 50;
            qp_scale_factor                         = 3;

            break;
        case 10:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 32;
            disallow_below_32x32_cost_th_multiplier = 2;
            disallow_below_64x64_cost_th_multiplier = 2;
            dev_16x16_to_8x8_th                     = 200;
            dev_32x32_to_16x16_th                   = 75;
            qp_scale_factor                         = 3;
            break;
        case 11:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 32;
            disallow_below_32x32_cost_th_multiplier = 2;
            disallow_below_64x64_cost_th_multiplier = 2;
            dev_16x16_to_8x8_th                     = 250;
            dev_32x32_to_16x16_th                   = 125;
            qp_scale_factor                         = 3;
            break;
        case 12:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 32;
            disallow_below_32x32_cost_th_multiplier = 4;
            disallow_below_64x64_cost_th_multiplier = 2;
            dev_16x16_to_8x8_th                     = 250;
            dev_32x32_to_16x16_th                   = 150;
            qp_scale_factor                         = 4;
            break;
        case 13:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 64;
            disallow_below_32x32_cost_th_multiplier = 4;
            disallow_below_64x64_cost_th_multiplier = 2;
            dev_16x16_to_8x8_th                     = 250;
            dev_32x32_to_16x16_th                   = 150;
            qp_scale_factor                         = 4;
            break;
        case 14:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 64;
            disallow_below_32x32_cost_th_multiplier = 4;
            disallow_below_64x64_cost_th_multiplier = 4;
            dev_16x16_to_8x8_th                     = 250;
            dev_32x32_to_16x16_th                   = 150;
            qp_scale_factor                         = 4;
            break;
        case 15:
            depth_removal_ctrls->enabled            = 1;
            disallow_below_16x16_cost_th_multiplier = 96;
            disallow_below_32x32_cost_th_multiplier = 6;
            disallow_below_64x64_cost_th_multiplier = 6;
            dev_16x16_to_8x8_th                     = 300;
            dev_32x32_to_16x16_th                   = 200;
            qp_scale_factor                         = 4;
            break;
        }
        if (depth_removal_ctrls->enabled) {
            //dev_16x16_to_8x8_th , dev_32x32_to_16x16_th = f(me_8x8_cost_variance)
            me_8x8_cost_variance /= MAX((MAX(63 - (pcs->picture_qp + 10), 1)), 1);
            if (me_8x8_cost_variance < LOW_8x8_DIST_VAR_TH) {
                dev_16x16_to_8x8_th = dev_16x16_to_8x8_th << 2;
            } else if (me_8x8_cost_variance < HIGH_8x8_DIST_VAR_TH) {
                dev_16x16_to_8x8_th   = dev_16x16_to_8x8_th << 1;
                dev_32x32_to_16x16_th = dev_32x32_to_16x16_th >> 1;
            } else {
                dev_16x16_to_8x8_th   = 0;
                dev_32x32_to_16x16_th = 0;
            }

            //dev_16x16_to_8x8_th , dev_32x32_to_16x16_th = f(QP)
            dev_16x16_to_8x8_th *= MAX((MAX(63 - (pcs->picture_qp + 10), 1) >> 4), 1) * qp_scale_factor;
            dev_32x32_to_16x16_th *= MAX((MAX(63 - (pcs->picture_qp + 10), 1) >> 4), 1) * qp_scale_factor;

            // dev_32x32_to_8x8_th = f(dev_32x32_to_16x16_th); a bit higher
            dev_32x32_to_8x8_th = (dev_32x32_to_16x16_th * ((1 << 2) + 1)) >> 2;

            uint64_t disallow_below_16x16_cost_th = disallow_below_16x16_cost_th_multiplier
                ? RDCOST(fast_lambda, cost_th_rate, (sb_size >> 1) * disallow_below_16x16_cost_th_multiplier)
                : 0;
            uint64_t disallow_below_32x32_cost_th = disallow_below_32x32_cost_th_multiplier
                ? RDCOST(fast_lambda, cost_th_rate, (sb_size >> 1) * disallow_below_32x32_cost_th_multiplier)
                : 0;
            uint64_t disallow_below_64x64_cost_th = disallow_below_64x64_cost_th_multiplier
                ? RDCOST(fast_lambda, cost_th_rate, (sb_size >> 1) * disallow_below_64x64_cost_th_multiplier)
                : 0;

            uint64_t cost_64x64 = RDCOST(fast_lambda, 0, pcs->ppcs->me_64x64_distortion[ctx->sb_index]);
            uint64_t cost_32x32 = RDCOST(fast_lambda, 0, pcs->ppcs->me_32x32_distortion[ctx->sb_index]);
            uint64_t cost_16x16 = RDCOST(fast_lambda, 0, pcs->ppcs->me_16x16_distortion[ctx->sb_index]);
            uint64_t cost_8x8   = RDCOST(fast_lambda, 0, pcs->ppcs->me_8x8_distortion[ctx->sb_index]);

            int64_t dev_32x32_to_16x16 = (int64_t)(((int64_t)MAX(cost_32x32, 1) - (int64_t)MAX(cost_16x16, 1)) * 1000) /
                (int64_t)MAX(cost_16x16, 1);

            int64_t dev_32x32_to_8x8 = (int64_t)(((int64_t)MAX(cost_32x32, 1) - (int64_t)MAX(cost_8x8, 1)) * 1000) /
                (int64_t)MAX(cost_8x8, 1);

            int64_t dev_16x16_to_8x8 = (int64_t)(((int64_t)MAX(cost_16x16, 1) - (int64_t)MAX(cost_8x8, 1)) * 1000) /
                (int64_t)MAX(cost_8x8, 1);
            depth_removal_ctrls->disallow_below_64x64 = (b64_geom->width % 64 == 0 && b64_geom->height % 64 == 0)
                ? (depth_removal_ctrls->disallow_below_64x64 || cost_64x64 < disallow_below_64x64_cost_th)
                : 0;

            depth_removal_ctrls->disallow_below_32x32 = (b64_geom->width % 32 == 0 && b64_geom->height % 32 == 0)
                ? (depth_removal_ctrls->disallow_below_32x32 || cost_32x32 < disallow_below_32x32_cost_th ||
                   (dev_32x32_to_16x16 < dev_32x32_to_16x16_th && dev_32x32_to_8x8 < dev_32x32_to_8x8_th))
                : 0;

            depth_removal_ctrls->disallow_below_16x16 = (b64_geom->width % 16 == 0 && b64_geom->height % 16 == 0)
                ? (depth_removal_ctrls->disallow_below_16x16 || cost_16x16 < disallow_below_16x16_cost_th ||
                   dev_16x16_to_8x8 < dev_16x16_to_8x8_th)
                : 0;
            if (!ctx->disallow_4x4) {
                uint64_t disallow_4x4_cost_th_multiplier = 64;
                uint64_t disallow_4x4_cost_th            = RDCOST(
                    fast_lambda, cost_th_rate, (sb_size >> 1) * disallow_4x4_cost_th_multiplier);

                if (cost_8x8 < disallow_4x4_cost_th && me_8x8_cost_variance < LOW_8x8_DIST_VAR_TH)
                    ctx->disallow_4x4 = 1;
            }
        }
    }
}
/*
 * Control NSQ search
 */
static void md_nsq_motion_search_controls(ModeDecisionContext *ctx, uint8_t md_nsq_mv_search_level) {
    MdNsqMotionSearchCtrls *md_nsq_me_ctrls = &ctx->md_nsq_me_ctrls;

    switch (md_nsq_mv_search_level) {
    case 0: md_nsq_me_ctrls->enabled = 0; break;
    case 1:
        md_nsq_me_ctrls->enabled                = 1;
        md_nsq_me_ctrls->use_ssd                = 0;
        md_nsq_me_ctrls->full_pel_search_width  = 31;
        md_nsq_me_ctrls->full_pel_search_height = 31;
        md_nsq_me_ctrls->enable_psad            = 1;
        break;

    case 2:
        md_nsq_me_ctrls->enabled                = 1;
        md_nsq_me_ctrls->use_ssd                = 0;
        md_nsq_me_ctrls->full_pel_search_width  = 15;
        md_nsq_me_ctrls->full_pel_search_height = 15;
        md_nsq_me_ctrls->enable_psad            = 1;
        break;
    case 3:
        md_nsq_me_ctrls->enabled                = 1;
        md_nsq_me_ctrls->use_ssd                = 0;
        md_nsq_me_ctrls->full_pel_search_width  = 11;
        md_nsq_me_ctrls->full_pel_search_height = 11;
        md_nsq_me_ctrls->enable_psad            = 1;
        break;
    case 4:
        md_nsq_me_ctrls->enabled                = 1;
        md_nsq_me_ctrls->use_ssd                = 0;
        md_nsq_me_ctrls->full_pel_search_width  = 8;
        md_nsq_me_ctrls->full_pel_search_height = 7;
        md_nsq_me_ctrls->enable_psad            = 1;
        break;
    default: assert(0); break;
    }
}
void svt_aom_md_pme_search_controls(ModeDecisionContext *ctx, uint8_t md_pme_level) {
    MdPmeCtrls *md_pme_ctrls = &ctx->md_pme_ctrls;

    switch (md_pme_level) {
    case 0: md_pme_ctrls->enabled = 0; break;
    case 1:
        md_pme_ctrls->enabled                       = 1;
        md_pme_ctrls->use_ssd                       = 1;
        md_pme_ctrls->full_pel_search_width         = 31;
        md_pme_ctrls->full_pel_search_height        = 31;
        md_pme_ctrls->early_check_mv_th_multiplier  = MIN_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th      = MAX_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th        = MIN_SIGNED_VALUE;
        md_pme_ctrls->post_fp_pme_to_me_cost_th     = MAX_SIGNED_VALUE;
        md_pme_ctrls->post_fp_pme_to_me_mv_th       = MIN_SIGNED_VALUE;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 0;
        md_pme_ctrls->enable_psad                   = 0;
        break;
    case 2:
        md_pme_ctrls->enabled                       = 1;
        md_pme_ctrls->use_ssd                       = 1;
        md_pme_ctrls->full_pel_search_width         = 7;
        md_pme_ctrls->full_pel_search_height        = 5;
        md_pme_ctrls->early_check_mv_th_multiplier  = MIN_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th      = MAX_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th        = MIN_SIGNED_VALUE;
        md_pme_ctrls->post_fp_pme_to_me_cost_th     = MAX_SIGNED_VALUE;
        md_pme_ctrls->post_fp_pme_to_me_mv_th       = MIN_SIGNED_VALUE;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 0;
        md_pme_ctrls->enable_psad                   = 0;
        break;
    case 3:
        md_pme_ctrls->enabled                       = 1;
        md_pme_ctrls->use_ssd                       = 1;
        md_pme_ctrls->full_pel_search_width         = 7;
        md_pme_ctrls->full_pel_search_height        = 5;
        md_pme_ctrls->early_check_mv_th_multiplier  = MIN_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th      = MAX_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th        = 4;
        md_pme_ctrls->post_fp_pme_to_me_cost_th     = MAX_SIGNED_VALUE;
        md_pme_ctrls->post_fp_pme_to_me_mv_th       = MIN_SIGNED_VALUE;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 0;
        md_pme_ctrls->enable_psad                   = 0;
        break;
    case 4:
        md_pme_ctrls->enabled                       = 1;
        md_pme_ctrls->use_ssd                       = 0;
        md_pme_ctrls->full_pel_search_width         = 3;
        md_pme_ctrls->full_pel_search_height        = 3;
        md_pme_ctrls->early_check_mv_th_multiplier  = MIN_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th      = MAX_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th        = 16;
        md_pme_ctrls->post_fp_pme_to_me_cost_th     = MAX_SIGNED_VALUE;
        md_pme_ctrls->post_fp_pme_to_me_mv_th       = 32;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 0;
        md_pme_ctrls->enable_psad                   = 1;
        break;
    case 5:
        md_pme_ctrls->enabled                       = 1;
        md_pme_ctrls->use_ssd                       = 0;
        md_pme_ctrls->full_pel_search_width         = 7;
        md_pme_ctrls->full_pel_search_height        = 5;
        md_pme_ctrls->early_check_mv_th_multiplier  = MIN_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th      = MAX_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th        = 16;
        md_pme_ctrls->post_fp_pme_to_me_cost_th     = 50;
        md_pme_ctrls->post_fp_pme_to_me_mv_th       = MIN_SIGNED_VALUE;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 0;
        md_pme_ctrls->enable_psad                   = 1;
        break;
    case 6:
        md_pme_ctrls->enabled                       = 1;
        md_pme_ctrls->use_ssd                       = 0;
        md_pme_ctrls->full_pel_search_width         = 7;
        md_pme_ctrls->full_pel_search_height        = 5;
        md_pme_ctrls->early_check_mv_th_multiplier  = MIN_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th      = 100;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th        = 16;
        md_pme_ctrls->post_fp_pme_to_me_cost_th     = 25;
        md_pme_ctrls->post_fp_pme_to_me_mv_th       = 32;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 0;
        md_pme_ctrls->enable_psad                   = 1;
        break;
    case 7:
        md_pme_ctrls->enabled                       = 1;
        md_pme_ctrls->use_ssd                       = 0;
        md_pme_ctrls->full_pel_search_width         = 7;
        md_pme_ctrls->full_pel_search_height        = 5;
        md_pme_ctrls->early_check_mv_th_multiplier  = MIN_SIGNED_VALUE;
        md_pme_ctrls->pre_fp_pme_to_me_cost_th      = 25;
        md_pme_ctrls->pre_fp_pme_to_me_mv_th        = 16;
        md_pme_ctrls->post_fp_pme_to_me_cost_th     = 5;
        md_pme_ctrls->post_fp_pme_to_me_mv_th       = 32;
        md_pme_ctrls->modulate_pme_for_blk_size_res = 0;
        md_pme_ctrls->enable_psad                   = 1;
        break;
    default: assert(0); break;
    }
}
static void set_subres_controls(ModeDecisionContext *ctx, uint8_t subres_level) {
    SubresCtrls *subres_ctrls = &ctx->subres_ctrls;

    switch (subres_level) {
    case 0: subres_ctrls->step = 0; break;
    case 1: subres_ctrls->step = 1; break;
    case 2: subres_ctrls->step = 2; break;
    default: assert(0); break;
    }
    // Set the TH used to determine if subres is safe to use (based on ODD vs. EVEN rows' distortion)
    if (subres_ctrls->step == 0)
        subres_ctrls->odd_to_even_deviation_th = 0;
    else
        subres_ctrls->odd_to_even_deviation_th = 5;
}
static void set_pf_controls(ModeDecisionContext *ctx, uint8_t pf_level) {
    PfCtrls *pf_ctrls = &ctx->pf_ctrls;

    switch (pf_level) {
    case 0: pf_ctrls->pf_shape = ONLY_DC_SHAPE; break;
    case 1: pf_ctrls->pf_shape = DEFAULT_SHAPE; break;
    case 2: pf_ctrls->pf_shape = N2_SHAPE; break;
    case 3: pf_ctrls->pf_shape = N4_SHAPE; break;
    default: assert(0); break;
    }
}
/*
 * Control Adaptive ME search
 */
static void md_sq_motion_search_controls(ModeDecisionContext *ctx, uint8_t md_sq_mv_search_level) {
    MdSqMotionSearchCtrls *md_sq_me_ctrls = &ctx->md_sq_me_ctrls;

    switch (md_sq_mv_search_level) {
    case 0: md_sq_me_ctrls->enabled = 0; break;
    case 1:
        md_sq_me_ctrls->enabled = 1;
        md_sq_me_ctrls->use_ssd = 0;

        md_sq_me_ctrls->pame_distortion_th = 10;

        md_sq_me_ctrls->sprs_lev0_enabled    = 1;
        md_sq_me_ctrls->sprs_lev0_step       = 4;
        md_sq_me_ctrls->sprs_lev0_w          = 15;
        md_sq_me_ctrls->sprs_lev0_h          = 15;
        md_sq_me_ctrls->max_sprs_lev0_w      = 150;
        md_sq_me_ctrls->max_sprs_lev0_h      = 150;
        md_sq_me_ctrls->sprs_lev0_multiplier = 500;

        md_sq_me_ctrls->sprs_lev1_enabled    = 1;
        md_sq_me_ctrls->sprs_lev1_step       = 2;
        md_sq_me_ctrls->sprs_lev1_w          = 4;
        md_sq_me_ctrls->sprs_lev1_h          = 4;
        md_sq_me_ctrls->max_sprs_lev1_w      = 50;
        md_sq_me_ctrls->max_sprs_lev1_h      = 50;
        md_sq_me_ctrls->sprs_lev1_multiplier = 500;

        md_sq_me_ctrls->sprs_lev2_enabled = 1;
        md_sq_me_ctrls->sprs_lev2_step    = 1;
        md_sq_me_ctrls->sprs_lev2_w       = 3;
        md_sq_me_ctrls->sprs_lev2_h       = 3;
        md_sq_me_ctrls->enable_psad       = 1;
        break;
    case 2:
        md_sq_me_ctrls->enabled = 1;
        md_sq_me_ctrls->use_ssd = 0;

        md_sq_me_ctrls->pame_distortion_th = 10;

        md_sq_me_ctrls->sprs_lev0_enabled    = 1;
        md_sq_me_ctrls->sprs_lev0_step       = 4;
        md_sq_me_ctrls->sprs_lev0_w          = 15;
        md_sq_me_ctrls->sprs_lev0_h          = 15;
        md_sq_me_ctrls->max_sprs_lev0_w      = 150;
        md_sq_me_ctrls->max_sprs_lev0_h      = 150;
        md_sq_me_ctrls->sprs_lev0_multiplier = 400;

        md_sq_me_ctrls->sprs_lev1_enabled    = 1;
        md_sq_me_ctrls->sprs_lev1_step       = 2;
        md_sq_me_ctrls->sprs_lev1_w          = 4;
        md_sq_me_ctrls->sprs_lev1_h          = 4;
        md_sq_me_ctrls->max_sprs_lev1_w      = 50;
        md_sq_me_ctrls->max_sprs_lev1_h      = 50;
        md_sq_me_ctrls->sprs_lev1_multiplier = 400;

        md_sq_me_ctrls->sprs_lev2_enabled = 1;
        md_sq_me_ctrls->sprs_lev2_step    = 1;
        md_sq_me_ctrls->sprs_lev2_w       = 3;
        md_sq_me_ctrls->sprs_lev2_h       = 3;
        md_sq_me_ctrls->enable_psad       = 1;
        break;
    case 3:
        md_sq_me_ctrls->enabled = 1;
        md_sq_me_ctrls->use_ssd = 0;

        md_sq_me_ctrls->pame_distortion_th = 10;

        md_sq_me_ctrls->sprs_lev0_enabled    = 1;
        md_sq_me_ctrls->sprs_lev0_step       = 4;
        md_sq_me_ctrls->sprs_lev0_w          = 15;
        md_sq_me_ctrls->sprs_lev0_h          = 15;
        md_sq_me_ctrls->max_sprs_lev0_w      = 150;
        md_sq_me_ctrls->max_sprs_lev0_h      = 150;
        md_sq_me_ctrls->sprs_lev0_multiplier = 300;

        md_sq_me_ctrls->sprs_lev1_enabled    = 1;
        md_sq_me_ctrls->sprs_lev1_step       = 2;
        md_sq_me_ctrls->sprs_lev1_w          = 4;
        md_sq_me_ctrls->sprs_lev1_h          = 4;
        md_sq_me_ctrls->max_sprs_lev1_w      = 50;
        md_sq_me_ctrls->max_sprs_lev1_h      = 50;
        md_sq_me_ctrls->sprs_lev1_multiplier = 300;

        md_sq_me_ctrls->sprs_lev2_enabled = 1;
        md_sq_me_ctrls->sprs_lev2_step    = 1;
        md_sq_me_ctrls->sprs_lev2_w       = 3;
        md_sq_me_ctrls->sprs_lev2_h       = 3;
        md_sq_me_ctrls->enable_psad       = 1;
        break;
    case 4:
        md_sq_me_ctrls->enabled            = 1;
        md_sq_me_ctrls->use_ssd            = 0;
        md_sq_me_ctrls->pame_distortion_th = 10;

        md_sq_me_ctrls->sprs_lev0_enabled    = 1;
        md_sq_me_ctrls->sprs_lev0_step       = 4;
        md_sq_me_ctrls->sprs_lev0_w          = 15;
        md_sq_me_ctrls->sprs_lev0_h          = 15;
        md_sq_me_ctrls->max_sprs_lev0_w      = 150;
        md_sq_me_ctrls->max_sprs_lev0_h      = 150;
        md_sq_me_ctrls->sprs_lev0_multiplier = 100;

        md_sq_me_ctrls->sprs_lev1_enabled    = 1;
        md_sq_me_ctrls->sprs_lev1_step       = 2;
        md_sq_me_ctrls->sprs_lev1_w          = 4;
        md_sq_me_ctrls->sprs_lev1_h          = 4;
        md_sq_me_ctrls->max_sprs_lev1_w      = 50;
        md_sq_me_ctrls->max_sprs_lev1_h      = 50;
        md_sq_me_ctrls->sprs_lev1_multiplier = 100;

        md_sq_me_ctrls->sprs_lev2_enabled = 1;
        md_sq_me_ctrls->sprs_lev2_step    = 1;
        md_sq_me_ctrls->sprs_lev2_w       = 3;
        md_sq_me_ctrls->sprs_lev2_h       = 3;
        md_sq_me_ctrls->enable_psad       = 1;
        break;
    default: assert(0); break;
    }
}
/*
 * Control Subpel search of ME MV(s)
 */
static void md_subpel_me_controls(ModeDecisionContext *ctx, uint8_t md_subpel_me_level, bool rtc_tune) {
    MdSubPelSearchCtrls *md_subpel_me_ctrls = &ctx->md_subpel_me_ctrls;

    switch (md_subpel_me_level) {
    case 0: md_subpel_me_ctrls->enabled = 0; break;
    case 1:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_8_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 2;
        md_subpel_me_ctrls->max_precision         = EIGHTH_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE;
        md_subpel_me_ctrls->pred_variance_th      = 0;
        md_subpel_me_ctrls->abs_th_mult           = 0;
        md_subpel_me_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement  = 0;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        md_subpel_me_ctrls->min_blk_sz            = 0;
        md_subpel_me_ctrls->mvp_th                = 0;
        break;
    case 2:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 2;
        md_subpel_me_ctrls->max_precision         = QUARTER_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE;
        md_subpel_me_ctrls->pred_variance_th      = 0;
        md_subpel_me_ctrls->abs_th_mult           = 0;
        md_subpel_me_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement  = 0;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        md_subpel_me_ctrls->min_blk_sz            = 4;
        md_subpel_me_ctrls->mvp_th                = 18;
        break;
    case 3:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 2;
        md_subpel_me_ctrls->max_precision         = QUARTER_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE;
        md_subpel_me_ctrls->pred_variance_th      = 0;
        md_subpel_me_ctrls->abs_th_mult           = 0;
        md_subpel_me_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement  = 0;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        md_subpel_me_ctrls->min_blk_sz            = 4;
        md_subpel_me_ctrls->mvp_th                = 12;
        break;
    case 4:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->max_precision         = QUARTER_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE;
        md_subpel_me_ctrls->pred_variance_th      = 0;
        md_subpel_me_ctrls->abs_th_mult           = 0;
        md_subpel_me_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement  = 0;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        md_subpel_me_ctrls->min_blk_sz            = 4;
        md_subpel_me_ctrls->mvp_th                = 12;
        break;
    case 5:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_8_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 2;
        md_subpel_me_ctrls->max_precision         = EIGHTH_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th      = 0;
        md_subpel_me_ctrls->abs_th_mult           = 0;
        md_subpel_me_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement  = 1;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        md_subpel_me_ctrls->min_blk_sz            = 4;
        md_subpel_me_ctrls->mvp_th                = 12;
        break;
    case 6:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_8_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 2;
        md_subpel_me_ctrls->max_precision         = EIGHTH_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th      = 0;
        md_subpel_me_ctrls->abs_th_mult           = 0;
        md_subpel_me_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement  = 2;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        md_subpel_me_ctrls->min_blk_sz            = 4;
        md_subpel_me_ctrls->mvp_th                = 12;
        break;
    case 7:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 2;
        md_subpel_me_ctrls->max_precision         = QUARTER_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th      = 0;
        md_subpel_me_ctrls->abs_th_mult           = 0;
        md_subpel_me_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement  = 2;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        md_subpel_me_ctrls->min_blk_sz            = 4;
        md_subpel_me_ctrls->mvp_th                = 12;
        break;
    case 8:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->max_precision         = QUARTER_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th      = 0;
        md_subpel_me_ctrls->abs_th_mult           = 0;
        md_subpel_me_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement  = 1;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        md_subpel_me_ctrls->min_blk_sz            = 4;
        md_subpel_me_ctrls->mvp_th                = 12;
        break;
    case 9:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->max_precision         = QUARTER_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th      = 0;
        md_subpel_me_ctrls->abs_th_mult           = 0;
        md_subpel_me_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement  = 2;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        md_subpel_me_ctrls->min_blk_sz            = 4;
        md_subpel_me_ctrls->mvp_th                = 12;
        break;
    case 10:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->max_precision         = QUARTER_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th      = 0;
        md_subpel_me_ctrls->abs_th_mult           = 0;
        md_subpel_me_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement  = 3;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        md_subpel_me_ctrls->min_blk_sz            = 4;
        md_subpel_me_ctrls->mvp_th                = 12;
        break;
    case 11:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->max_precision         = QUARTER_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th      = 0;
        md_subpel_me_ctrls->abs_th_mult           = 0;
        md_subpel_me_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement  = 2;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        md_subpel_me_ctrls->min_blk_sz            = 4;
        md_subpel_me_ctrls->mvp_th                = 12;
        break;
    case 12:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->max_precision         = QUARTER_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th      = 0;
        md_subpel_me_ctrls->abs_th_mult           = 0;
        md_subpel_me_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement  = 3;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        md_subpel_me_ctrls->min_blk_sz            = 4;
        md_subpel_me_ctrls->mvp_th                = 12;
        break;
    case 13:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->max_precision         = QUARTER_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th      = 50;
        md_subpel_me_ctrls->abs_th_mult           = 2;
        md_subpel_me_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_me_ctrls->skip_diag_refinement  = 3;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        md_subpel_me_ctrls->min_blk_sz            = 4;
        md_subpel_me_ctrls->mvp_th                = 12;
        break;
    case 14:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->max_precision         = QUARTER_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th      = 100;
        md_subpel_me_ctrls->abs_th_mult           = 2;
        md_subpel_me_ctrls->round_dev_th          = rtc_tune ? MAX_SIGNED_VALUE : -25;
        md_subpel_me_ctrls->skip_diag_refinement  = 4;
        md_subpel_me_ctrls->skip_zz_mv            = 0;
        md_subpel_me_ctrls->min_blk_sz            = 4;
        md_subpel_me_ctrls->mvp_th                = 12;
        break;
    case 15:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->max_precision         = QUARTER_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th      = 100;
        md_subpel_me_ctrls->abs_th_mult           = 2;
        md_subpel_me_ctrls->round_dev_th          = rtc_tune ? MAX_SIGNED_VALUE : -25;
        md_subpel_me_ctrls->skip_diag_refinement  = 4;
        md_subpel_me_ctrls->skip_zz_mv            = rtc_tune ? 0 : 1;
        md_subpel_me_ctrls->min_blk_sz            = 4;
        md_subpel_me_ctrls->mvp_th                = 12;
        break;
    case 16:
        md_subpel_me_ctrls->enabled               = 1;
        md_subpel_me_ctrls->subpel_search_type    = USE_4_TAPS;
        md_subpel_me_ctrls->subpel_iters_per_step = 1;
        md_subpel_me_ctrls->max_precision         = HALF_PEL;
        md_subpel_me_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_me_ctrls->pred_variance_th      = 100;
        md_subpel_me_ctrls->abs_th_mult           = 2;
        md_subpel_me_ctrls->round_dev_th          = -25;
        md_subpel_me_ctrls->skip_diag_refinement  = 4;
        md_subpel_me_ctrls->skip_zz_mv            = 1;
        md_subpel_me_ctrls->min_blk_sz            = 4;
        md_subpel_me_ctrls->mvp_th                = 12;
        break;

    default: assert(0); break;
    }
}

/*
 * Control Subpel search of PME MV(s)
 */
static void md_subpel_pme_controls(ModeDecisionContext *ctx, uint8_t md_subpel_pme_level) {
    MdSubPelSearchCtrls *md_subpel_pme_ctrls = &ctx->md_subpel_pme_ctrls;

    switch (md_subpel_pme_level) {
    case 0: md_subpel_pme_ctrls->enabled = 0; break;
    case 1:
        md_subpel_pme_ctrls->enabled               = 1;
        md_subpel_pme_ctrls->subpel_search_type    = USE_8_TAPS;
        md_subpel_pme_ctrls->subpel_iters_per_step = 2;
        md_subpel_pme_ctrls->max_precision         = EIGHTH_PEL;
        md_subpel_pme_ctrls->subpel_search_method  = SUBPEL_TREE;
        md_subpel_pme_ctrls->pred_variance_th      = 0;
        md_subpel_pme_ctrls->abs_th_mult           = 0;
        md_subpel_pme_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_pme_ctrls->skip_zz_mv            = 0;
        md_subpel_pme_ctrls->min_blk_sz            = 0;
        md_subpel_pme_ctrls->mvp_th                = 0;
        break;
    case 2:
        md_subpel_pme_ctrls->enabled               = 1;
        md_subpel_pme_ctrls->subpel_search_type    = USE_8_TAPS;
        md_subpel_pme_ctrls->subpel_iters_per_step = 2;
        md_subpel_pme_ctrls->max_precision         = EIGHTH_PEL;
        md_subpel_pme_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_pme_ctrls->pred_variance_th      = 0;
        md_subpel_pme_ctrls->abs_th_mult           = 0;
        md_subpel_pme_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_pme_ctrls->skip_zz_mv            = 0;
        md_subpel_pme_ctrls->min_blk_sz            = 0;
        md_subpel_pme_ctrls->mvp_th                = 0;
        break;
    case 3:
        md_subpel_pme_ctrls->enabled               = 1;
        md_subpel_pme_ctrls->subpel_search_type    = USE_8_TAPS;
        md_subpel_pme_ctrls->subpel_iters_per_step = 2;
        md_subpel_pme_ctrls->max_precision         = EIGHTH_PEL;
        md_subpel_pme_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_pme_ctrls->pred_variance_th      = 0;
        md_subpel_pme_ctrls->abs_th_mult           = 0;
        md_subpel_pme_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_pme_ctrls->skip_zz_mv            = 0;
        md_subpel_pme_ctrls->min_blk_sz            = 0;
        md_subpel_pme_ctrls->mvp_th                = 0;
        break;
    case 4:
        md_subpel_pme_ctrls->enabled               = 1;
        md_subpel_pme_ctrls->subpel_search_type    = USE_8_TAPS;
        md_subpel_pme_ctrls->subpel_iters_per_step = 2;
        md_subpel_pme_ctrls->max_precision         = HALF_PEL;
        md_subpel_pme_ctrls->subpel_search_method  = SUBPEL_TREE_PRUNED;
        md_subpel_pme_ctrls->pred_variance_th      = 0;
        md_subpel_pme_ctrls->abs_th_mult           = 0;
        md_subpel_pme_ctrls->round_dev_th          = MAX_SIGNED_VALUE;
        md_subpel_pme_ctrls->skip_zz_mv            = 0;
        md_subpel_pme_ctrls->min_blk_sz            = 0;
        md_subpel_pme_ctrls->mvp_th                = 0;
        break;
    default: assert(0); break;
    }
}
/*
 * Control RDOQ
 */
static void set_rdoq_controls(ModeDecisionContext *ctx, uint8_t rdoq_level) {
    RdoqCtrls *rdoq_ctrls = &ctx->rdoq_ctrls;

    switch (rdoq_level) {
    case 0: rdoq_ctrls->enabled = 0; break;
    case 1:
        rdoq_ctrls->enabled          = 1;
        rdoq_ctrls->eob_fast_y_inter = 0;
        rdoq_ctrls->eob_fast_y_intra = 0;
#if TUNE_CHROMA_SSIM
        rdoq_ctrls->eob_fast_uv_inter = 0;
#else
        rdoq_ctrls->eob_fast_uv_inter = 1;
#endif
        rdoq_ctrls->eob_fast_uv_intra = 0;
        rdoq_ctrls->fp_q_y            = 1;
        rdoq_ctrls->fp_q_uv           = 1;
        rdoq_ctrls->satd_factor       = (uint8_t)~0;
        rdoq_ctrls->early_exit_th     = 0;
        rdoq_ctrls->skip_uv           = 0;
        rdoq_ctrls->dct_dct_only      = 0;
        rdoq_ctrls->eob_th            = (uint8_t)~0;
        rdoq_ctrls->eob_fast_th       = (uint8_t)~0;
        break;
    case 2:
        rdoq_ctrls->enabled          = 1;
        rdoq_ctrls->eob_fast_y_inter = 0;
        rdoq_ctrls->eob_fast_y_intra = 0;
#if TUNE_CHROMA_SSIM
        rdoq_ctrls->eob_fast_uv_inter = 0;
#else
        rdoq_ctrls->eob_fast_uv_inter = 1;
#endif
        rdoq_ctrls->eob_fast_uv_intra = 0;
        rdoq_ctrls->fp_q_y            = 1;
        rdoq_ctrls->fp_q_uv           = 1;
        rdoq_ctrls->satd_factor       = (uint8_t)~0;
        rdoq_ctrls->early_exit_th     = 0;
        rdoq_ctrls->skip_uv           = 1;
        rdoq_ctrls->dct_dct_only      = 0;
        rdoq_ctrls->eob_th            = (uint8_t)~0;
        rdoq_ctrls->eob_fast_th       = (uint8_t)~0;
        break;
    case 3:
        rdoq_ctrls->enabled          = 1;
        rdoq_ctrls->eob_fast_y_inter = 0;
        rdoq_ctrls->eob_fast_y_intra = 0;
#if TUNE_CHROMA_SSIM
        rdoq_ctrls->eob_fast_uv_inter = 0;
#else
        rdoq_ctrls->eob_fast_uv_inter = 1;
#endif
        rdoq_ctrls->eob_fast_uv_intra = 0;
        rdoq_ctrls->fp_q_y            = 1;
        rdoq_ctrls->fp_q_uv           = 1;
        rdoq_ctrls->satd_factor       = (uint8_t)~0;
        rdoq_ctrls->early_exit_th     = 0;
        rdoq_ctrls->skip_uv           = 1;
        rdoq_ctrls->dct_dct_only      = 1;
        rdoq_ctrls->eob_th            = (uint8_t)~0;
        rdoq_ctrls->eob_fast_th       = (uint8_t)~0;
        break;
    case 4:
        rdoq_ctrls->enabled          = 1;
        rdoq_ctrls->eob_fast_y_inter = 0;
        rdoq_ctrls->eob_fast_y_intra = 0;
#if TUNE_CHROMA_SSIM
        rdoq_ctrls->eob_fast_uv_inter = 0;
#else
        rdoq_ctrls->eob_fast_uv_inter = 1;
#endif
        rdoq_ctrls->eob_fast_uv_intra = 0;
        rdoq_ctrls->fp_q_y            = 1;
        rdoq_ctrls->fp_q_uv           = 1;
        rdoq_ctrls->satd_factor       = (uint8_t)~0;
        rdoq_ctrls->early_exit_th     = 0;
        rdoq_ctrls->skip_uv           = 1;
        rdoq_ctrls->dct_dct_only      = 1;
        rdoq_ctrls->eob_th            = (uint8_t)~0;
        rdoq_ctrls->eob_fast_th       = 30;
        break;
    case 5:
        rdoq_ctrls->enabled          = 1;
        rdoq_ctrls->eob_fast_y_inter = 0;
        rdoq_ctrls->eob_fast_y_intra = 0;
#if TUNE_CHROMA_SSIM
        rdoq_ctrls->eob_fast_uv_inter = 0;
#else
        rdoq_ctrls->eob_fast_uv_inter = 1;
#endif
        rdoq_ctrls->eob_fast_uv_intra = 0;
        rdoq_ctrls->fp_q_y            = 1;
        rdoq_ctrls->fp_q_uv           = 1;
        rdoq_ctrls->satd_factor       = (uint8_t)~0;
        rdoq_ctrls->early_exit_th     = 0;
        rdoq_ctrls->skip_uv           = 1;
        rdoq_ctrls->dct_dct_only      = 1;
        rdoq_ctrls->eob_th            = 85;
        rdoq_ctrls->eob_fast_th       = 0;
        break;
    default: assert(0); break;
    }
}
/*
 * Settings for the parent SQ coeff-area based cycles reduction algorithm.
 */
static void set_parent_sq_coeff_area_based_cycles_reduction_ctrls(ModeDecisionContext *ctx, uint8_t resolution,
                                                                  uint8_t cycles_alloc_lvl) {
    ParentSqCmplxCtrls *cycle_red_ctrls = &ctx->psq_cplx_ctrls;
    switch (cycles_alloc_lvl) {
    case 0: cycle_red_ctrls->enabled = 0; break;
    case 1:
        cycle_red_ctrls->enabled = 1;

        // High frequency band THs/actions
        cycle_red_ctrls->high_freq_band1_th    = resolution <= INPUT_SIZE_360p_RANGE ? UNUSED_HIGH_FREQ_BAND_TH : 90;
        cycle_red_ctrls->high_freq_band1_level = 3;
        cycle_red_ctrls->high_freq_band2_th    = resolution <= INPUT_SIZE_360p_RANGE ? UNUSED_HIGH_FREQ_BAND_TH : 70;
        cycle_red_ctrls->high_freq_band2_level = 2;
        cycle_red_ctrls->high_freq_band3_th    = UNUSED_HIGH_FREQ_BAND_TH;
        cycle_red_ctrls->high_freq_band3_level = 0;

        // Low frequency band THs/actions
        cycle_red_ctrls->enable_zero_coeff_action = 1;
        cycle_red_ctrls->zero_coeff_action        = 1;
        cycle_red_ctrls->enable_one_coeff_action  = 0;
        cycle_red_ctrls->one_coeff_action         = 0;

        break;
    case 2:
        cycle_red_ctrls->enabled = 1;

        // High frequency band THs/actions
        cycle_red_ctrls->high_freq_band1_th    = resolution <= INPUT_SIZE_360p_RANGE ? UNUSED_HIGH_FREQ_BAND_TH : 90;
        cycle_red_ctrls->high_freq_band1_level = 3;
        cycle_red_ctrls->high_freq_band2_th    = resolution <= INPUT_SIZE_360p_RANGE ? UNUSED_HIGH_FREQ_BAND_TH : 70;
        cycle_red_ctrls->high_freq_band2_level = 2;
        cycle_red_ctrls->high_freq_band3_th    = UNUSED_HIGH_FREQ_BAND_TH;
        cycle_red_ctrls->high_freq_band3_level = 0;

        // Low frequency band THs/actions
        cycle_red_ctrls->enable_zero_coeff_action = 1;
        cycle_red_ctrls->zero_coeff_action        = 2;
        cycle_red_ctrls->enable_one_coeff_action  = 0;
        cycle_red_ctrls->one_coeff_action         = 0;
        break;
    case 3:
        cycle_red_ctrls->enabled = 1;
        // High frequency band THs/actions
        cycle_red_ctrls->high_freq_band1_th    = 90;
        cycle_red_ctrls->high_freq_band1_level = 3;
        cycle_red_ctrls->high_freq_band2_th    = 70;
        cycle_red_ctrls->high_freq_band2_level = 2;
        cycle_red_ctrls->high_freq_band3_th    = 50;
        cycle_red_ctrls->high_freq_band3_level = 1;
        // Low frequency band THs/actions
        cycle_red_ctrls->enable_zero_coeff_action = 1;
        cycle_red_ctrls->zero_coeff_action        = 3;
        cycle_red_ctrls->enable_one_coeff_action  = 1;
        cycle_red_ctrls->one_coeff_action         = 1;
        break;

    case 4:

        cycle_red_ctrls->enabled = 1;
        // High frequency band THs/actions
        cycle_red_ctrls->high_freq_band1_th    = 90;
        cycle_red_ctrls->high_freq_band1_level = 0;
        cycle_red_ctrls->high_freq_band2_th    = 70;
        cycle_red_ctrls->high_freq_band2_level = 3;
        cycle_red_ctrls->high_freq_band3_th    = 50;
        cycle_red_ctrls->high_freq_band3_level = 2;

        // Low frequency band THs/actions
        cycle_red_ctrls->enable_zero_coeff_action = 1;
        cycle_red_ctrls->zero_coeff_action        = 3;
        cycle_red_ctrls->enable_one_coeff_action  = 1;
        cycle_red_ctrls->one_coeff_action         = 1;

        break;
    case 5:
        cycle_red_ctrls->enabled = 1;

        // High frequency band THs/actions
        cycle_red_ctrls->high_freq_band1_th    = 90;
        cycle_red_ctrls->high_freq_band1_level = 0;
        cycle_red_ctrls->high_freq_band2_th    = 70;
        cycle_red_ctrls->high_freq_band2_level = 3;
        cycle_red_ctrls->high_freq_band3_th    = 50;
        cycle_red_ctrls->high_freq_band3_level = 2;

        // Low frequency band THs/actions
        cycle_red_ctrls->enable_zero_coeff_action = 1;
        cycle_red_ctrls->zero_coeff_action        = 0;
        cycle_red_ctrls->enable_one_coeff_action  = 1;
        cycle_red_ctrls->one_coeff_action         = 1;
        break;
    default: assert(0); break;
    }
}
static void set_sq_txs_ctrls(ModeDecisionContext *ctx, uint8_t psq_txs_lvl) {
    NsqPsqTxsCtrls *nsq_psq_txs_ctrls = &ctx->nsq_psq_txs_ctrls;
    switch (psq_txs_lvl) {
    case 0: nsq_psq_txs_ctrls->enabled = 0; break;
    case 1:
        nsq_psq_txs_ctrls->enabled     = 1;
        nsq_psq_txs_ctrls->hv_to_sq_th = 250;
        nsq_psq_txs_ctrls->h_to_v_th   = 25;
        break;
    case 2:
        nsq_psq_txs_ctrls->enabled     = 1;
        nsq_psq_txs_ctrls->hv_to_sq_th = 150;
        nsq_psq_txs_ctrls->h_to_v_th   = 15;
        break;
    default: assert(0); break;
    }
}

static void set_sq_pred_ctrls(ModeDecisionContext *ctx, uint8_t psq_pred_lvl) {
    NsqPsqPredCtrls *nsq_psq_pred_ctrls = &ctx->nsq_psq_pred_ctrls;
    switch (psq_pred_lvl) {
    case 0: nsq_psq_pred_ctrls->enabled = 0; break;
    case 1:
        nsq_psq_pred_ctrls->enabled = 1;
        nsq_psq_pred_ctrls->coef_th = 0;
        nsq_psq_pred_ctrls->cost_th = 20;
        break;
    case 2:
        nsq_psq_pred_ctrls->enabled = 1;
        nsq_psq_pred_ctrls->coef_th = 5;
        nsq_psq_pred_ctrls->cost_th = 10;
        break;
    case 3:
        nsq_psq_pred_ctrls->enabled = 1;
        nsq_psq_pred_ctrls->coef_th = 15;
        nsq_psq_pred_ctrls->cost_th = 10;
        break;
    case 4:
        nsq_psq_pred_ctrls->enabled = 1;
        nsq_psq_pred_ctrls->coef_th = 20;
        nsq_psq_pred_ctrls->cost_th = 10;
        break;
    default: assert(0); break;
    }
}

void svt_aom_set_txt_controls(ModeDecisionContext *ctx, uint8_t txt_level) {
    TxtControls *txt_ctrls = &ctx->txt_ctrls;

    switch (txt_level) {
    case 0:
        txt_ctrls->enabled = 0;

        txt_ctrls->txt_group_inter_lt_16x16    = 1;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 1;

        txt_ctrls->txt_group_intra_lt_16x16    = 1;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = 1;
        txt_ctrls->early_exit_dist_th          = 0;
        txt_ctrls->early_exit_coeff_th         = 0;
        txt_ctrls->satd_early_exit_th_intra    = 0;
        txt_ctrls->satd_early_exit_th_inter    = 0;
        txt_ctrls->txt_rate_cost_th            = 0;
        break;
    case 1:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16    = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = MAX_TX_TYPE_GROUP;

        txt_ctrls->txt_group_intra_lt_16x16    = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = MAX_TX_TYPE_GROUP;
        txt_ctrls->early_exit_dist_th          = 0;
        txt_ctrls->early_exit_coeff_th         = 0;
        txt_ctrls->satd_early_exit_th_intra    = 0;
        txt_ctrls->satd_early_exit_th_inter    = 0;
        txt_ctrls->txt_rate_cost_th            = 0;
        break;
    case 2:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16    = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_intra_lt_16x16    = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = MAX_TX_TYPE_GROUP;
        txt_ctrls->early_exit_dist_th          = 0;
        txt_ctrls->early_exit_coeff_th         = 0;
        txt_ctrls->satd_early_exit_th_intra    = 20;
        txt_ctrls->satd_early_exit_th_inter    = 20;
        txt_ctrls->txt_rate_cost_th            = 120;
        break;
    case 3:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16    = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = MAX_TX_TYPE_GROUP;

        txt_ctrls->txt_group_intra_lt_16x16    = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = MAX_TX_TYPE_GROUP;
        txt_ctrls->early_exit_dist_th          = 0;
        txt_ctrls->early_exit_coeff_th         = 0;
        txt_ctrls->satd_early_exit_th_intra    = 10;
        txt_ctrls->satd_early_exit_th_inter    = 10;
        txt_ctrls->txt_rate_cost_th            = 120;
        break;
    case 4:
        txt_ctrls->enabled                     = 1;
        txt_ctrls->txt_group_inter_lt_16x16    = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = MAX_TX_TYPE_GROUP;

        txt_ctrls->txt_group_intra_lt_16x16    = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = MAX_TX_TYPE_GROUP;
        txt_ctrls->early_exit_dist_th          = 0;
        txt_ctrls->early_exit_coeff_th         = 0;
        txt_ctrls->satd_early_exit_th_intra    = 10;
        txt_ctrls->satd_early_exit_th_inter    = 5;
        txt_ctrls->txt_rate_cost_th            = 100;
        break;
    case 5:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16    = 5;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 5;

        txt_ctrls->txt_group_intra_lt_16x16    = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = MAX_TX_TYPE_GROUP;
        txt_ctrls->early_exit_dist_th          = 0;
        txt_ctrls->early_exit_coeff_th         = 0;
        txt_ctrls->satd_early_exit_th_intra    = 10;
        txt_ctrls->satd_early_exit_th_inter    = 5;
        txt_ctrls->txt_rate_cost_th            = 100;
        break;
    case 6:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16    = 5;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 3;

        txt_ctrls->txt_group_intra_lt_16x16    = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = MAX_TX_TYPE_GROUP;
        txt_ctrls->early_exit_dist_th          = 0;
        txt_ctrls->early_exit_coeff_th         = 0;
        txt_ctrls->satd_early_exit_th_intra    = 10;
        txt_ctrls->satd_early_exit_th_inter    = 5;
        txt_ctrls->txt_rate_cost_th            = 100;
        break;
    case 7:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16    = 3;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 2;

        txt_ctrls->txt_group_intra_lt_16x16    = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = 4;
        txt_ctrls->early_exit_dist_th          = 0;
        txt_ctrls->early_exit_coeff_th         = 0;
        txt_ctrls->satd_early_exit_th_intra    = 10;
        txt_ctrls->satd_early_exit_th_inter    = 5;
        txt_ctrls->txt_rate_cost_th            = 100;
        break;
    case 8:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16    = 3;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 2;

        txt_ctrls->txt_group_intra_lt_16x16    = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = 4;

        txt_ctrls->early_exit_dist_th       = 100;
        txt_ctrls->early_exit_coeff_th      = 4;
        txt_ctrls->satd_early_exit_th_intra = 10;
        txt_ctrls->satd_early_exit_th_inter = 5;
        txt_ctrls->txt_rate_cost_th         = 100;
        break;
    case 9:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16    = 3;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 2;

        txt_ctrls->txt_group_intra_lt_16x16    = MAX_TX_TYPE_GROUP;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = 3;

        txt_ctrls->early_exit_dist_th       = 400;
        txt_ctrls->early_exit_coeff_th      = 8;
        txt_ctrls->satd_early_exit_th_intra = 10;
        txt_ctrls->satd_early_exit_th_inter = 5;
        txt_ctrls->txt_rate_cost_th         = 100;
        break;
    case 10:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16    = 3;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 2;

        txt_ctrls->txt_group_intra_lt_16x16    = 4;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = 2;
        txt_ctrls->early_exit_dist_th          = 400;
        txt_ctrls->early_exit_coeff_th         = 8;
        txt_ctrls->satd_early_exit_th_intra    = 10;
        txt_ctrls->satd_early_exit_th_inter    = 5;
        txt_ctrls->txt_rate_cost_th            = 100;
        break;
    case 11:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16    = 3;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 1;

        txt_ctrls->txt_group_intra_lt_16x16    = 4;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = 2;
        txt_ctrls->early_exit_dist_th          = 400;
        txt_ctrls->early_exit_coeff_th         = 8;
        txt_ctrls->satd_early_exit_th_intra    = 10;
        txt_ctrls->satd_early_exit_th_inter    = 5;
        txt_ctrls->txt_rate_cost_th            = 100;
        break;
    case 12:
        txt_ctrls->enabled = 1;

        txt_ctrls->txt_group_inter_lt_16x16    = 2;
        txt_ctrls->txt_group_inter_gt_eq_16x16 = 1;
        txt_ctrls->txt_group_intra_lt_16x16    = 3;
        txt_ctrls->txt_group_intra_gt_eq_16x16 = 2;
        txt_ctrls->early_exit_dist_th          = 400;
        txt_ctrls->early_exit_coeff_th         = 8;
        txt_ctrls->satd_early_exit_th_intra    = 10;
        txt_ctrls->satd_early_exit_th_inter    = 5;
        txt_ctrls->txt_rate_cost_th            = 100;
        break;
    default: assert(0); break;
    }
}
static void set_interpolation_search_level_ctrls(ModeDecisionContext *ctx, uint8_t interpolation_search_level) {
    InterpolationSearchCtrls *ifs_ctrls = &ctx->ifs_ctrls;

    switch (interpolation_search_level) {
    case 0:
        ifs_ctrls->level                 = IFS_OFF;
        ifs_ctrls->subsampled_distortion = 0;
        ifs_ctrls->skip_sse_rd_model     = 0;
        break;
    case 1:
        ifs_ctrls->level                 = IFS_MDS0;
        ifs_ctrls->subsampled_distortion = 0;
        ifs_ctrls->skip_sse_rd_model     = 0;
        break;
    case 2:
        ifs_ctrls->level                 = IFS_MDS1;
        ifs_ctrls->subsampled_distortion = 0;
        ifs_ctrls->skip_sse_rd_model     = 0;
        break;
    case 3:
        ifs_ctrls->level                 = IFS_MDS2;
        ifs_ctrls->subsampled_distortion = 0;
        ifs_ctrls->skip_sse_rd_model     = 0;
        break;
    case 4:
        ifs_ctrls->level                 = IFS_MDS3;
        ifs_ctrls->subsampled_distortion = 0;
        ifs_ctrls->skip_sse_rd_model     = 0;
        break;
    case 5:
        ifs_ctrls->level                 = IFS_MDS3;
        ifs_ctrls->subsampled_distortion = 1;
        ifs_ctrls->skip_sse_rd_model     = 1;
        break;
    default: assert(0); break;
    }
}
static void set_cand_reduction_ctrls(PictureControlSet *pcs, ModeDecisionContext *ctx, uint8_t cand_reduction_level,
                                     const uint32_t picture_qp, uint32_t me_8x8_cost_variance,
                                     uint32_t me_64x64_distortion, uint8_t l0_was_skip, uint8_t l1_was_skip,
                                     uint8_t ref_skip_perc) {
    CandReductionCtrls *cand_reduction_ctrls = &ctx->cand_reduction_ctrls;
    const bool          rtc_tune = (pcs->scs->static_config.pred_structure == SVT_AV1_PRED_LOW_DELAY_B) ? true : false;

    switch (cand_reduction_level) {
    case 0:

        // angular reduction at mds0
        cand_reduction_ctrls->mds0_reduce_intra = 0;

        // merge_inter_classes
        cand_reduction_ctrls->merge_inter_classes = 0;

        // redundant_cand_level
        cand_reduction_ctrls->redundant_cand_ctrls.score_th = 0;

        // use_neighbouring_mode
        cand_reduction_ctrls->use_neighbouring_mode_ctrls.enabled = 0;

        // near_count_ctrls
        cand_reduction_ctrls->near_count_ctrls.enabled         = 1;
        cand_reduction_ctrls->near_count_ctrls.near_count      = 3;
        cand_reduction_ctrls->near_count_ctrls.near_near_count = 3;

        // lpd1_mvp_best_me_list (LPD1 only signal)
        cand_reduction_ctrls->lpd1_mvp_best_me_list = 0;

        // cand_elimination_ctrls
        cand_reduction_ctrls->cand_elimination_ctrls.enabled = 0;

        // reduce_unipred_candidates
        cand_reduction_ctrls->reduce_unipred_candidates = 0;

        break;

    case 1:

        // angular reduction at mds0
        cand_reduction_ctrls->mds0_reduce_intra = 1;

        // merge_inter_classes
        cand_reduction_ctrls->merge_inter_classes = 0;

        // redundant_cand_level
        cand_reduction_ctrls->redundant_cand_ctrls.score_th = 0;

        // use_neighbouring_mode
        cand_reduction_ctrls->use_neighbouring_mode_ctrls.enabled = 0;

        // near_count_ctrls
        cand_reduction_ctrls->near_count_ctrls.enabled         = 1;
        cand_reduction_ctrls->near_count_ctrls.near_count      = 3;
        cand_reduction_ctrls->near_count_ctrls.near_near_count = 3;

        // lpd1_mvp_best_me_list (LPD1 only signal)
        cand_reduction_ctrls->lpd1_mvp_best_me_list = 0;

        // cand_elimination_ctrls
        cand_reduction_ctrls->cand_elimination_ctrls.enabled = 0;

        // reduce_unipred_candidates
        cand_reduction_ctrls->reduce_unipred_candidates = 0;

        break;
    case 2:

        // angular reduction at mds0
        cand_reduction_ctrls->mds0_reduce_intra = 1;

        // merge_inter_classes
        cand_reduction_ctrls->merge_inter_classes = 1;

        // redundant_cand_level
        cand_reduction_ctrls->redundant_cand_ctrls.score_th = 0;

        // use_neighbouring_mode
        cand_reduction_ctrls->use_neighbouring_mode_ctrls.enabled = 0;

        // near_count_ctrls
        cand_reduction_ctrls->near_count_ctrls.enabled         = 1;
        cand_reduction_ctrls->near_count_ctrls.near_count      = 3;
        cand_reduction_ctrls->near_count_ctrls.near_near_count = 3;

        // lpd1_mvp_best_me_list (LPD1 only signal)
        cand_reduction_ctrls->lpd1_mvp_best_me_list = 0;

        // cand_elimination_ctrls
        cand_reduction_ctrls->cand_elimination_ctrls.enabled         = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.dc_only         = 0;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_me   = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_pme  = 0;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_warp = 0;
        cand_reduction_ctrls->cand_elimination_ctrls.th_multiplier   = 1;

        // reduce_unipred_candidates
        cand_reduction_ctrls->reduce_unipred_candidates = 0;

        break;
    case 3:

        // angular reduction at mds0
        cand_reduction_ctrls->mds0_reduce_intra = 1;

        // merge_inter_classes
        cand_reduction_ctrls->merge_inter_classes = 1;

        // redundant_cand_level
        cand_reduction_ctrls->redundant_cand_ctrls.score_th = 0;

        // use_neighbouring_mode
        cand_reduction_ctrls->use_neighbouring_mode_ctrls.enabled = rtc_tune ? 0 : 1;

        // near_count_ctrls
        cand_reduction_ctrls->near_count_ctrls.enabled         = 1;
        cand_reduction_ctrls->near_count_ctrls.near_count      = 3;
        cand_reduction_ctrls->near_count_ctrls.near_near_count = 3;

        // lpd1_mvp_best_me_list (LPD1 only signal)
        cand_reduction_ctrls->lpd1_mvp_best_me_list = 0;

        // cand_elimination_ctrls
        cand_reduction_ctrls->cand_elimination_ctrls.enabled         = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.dc_only         = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_me   = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_pme  = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_warp = 0;
        cand_reduction_ctrls->cand_elimination_ctrls.th_multiplier   = 1;

        // reduce_unipred_candidates
        cand_reduction_ctrls->reduce_unipred_candidates = 0;

        break;

    case 4:

        // angular reduction at mds0
        cand_reduction_ctrls->mds0_reduce_intra = 1;

        // merge_inter_classes
        cand_reduction_ctrls->merge_inter_classes = 1;

        // redundant_cand_level
        cand_reduction_ctrls->redundant_cand_ctrls.score_th = 0;

        // use_neighbouring_mode
        cand_reduction_ctrls->use_neighbouring_mode_ctrls.enabled = rtc_tune ? 0 : 1;

        // near_count_ctrls
        cand_reduction_ctrls->near_count_ctrls.enabled         = 1;
        cand_reduction_ctrls->near_count_ctrls.near_count      = 1;
        cand_reduction_ctrls->near_count_ctrls.near_near_count = 3;

        // lpd1_mvp_best_me_list (LPD1 only signal)
        cand_reduction_ctrls->lpd1_mvp_best_me_list = 0;

        // cand_elimination_ctrls
        cand_reduction_ctrls->cand_elimination_ctrls.enabled         = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.dc_only         = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_me   = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_pme  = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_warp = 0;
        cand_reduction_ctrls->cand_elimination_ctrls.th_multiplier   = 1;

        // reduce_unipred_candidates
        cand_reduction_ctrls->reduce_unipred_candidates = 1;

        break;

    case 5:

        // angular reduction at mds0
        cand_reduction_ctrls->mds0_reduce_intra = 1;

        // merge_inter_classes
        cand_reduction_ctrls->merge_inter_classes = 1;

        // redundant_cand_level
        cand_reduction_ctrls->redundant_cand_ctrls.score_th = 8;
        cand_reduction_ctrls->redundant_cand_ctrls.mag_th   = 64;

        // use_neighbouring_mode
        cand_reduction_ctrls->use_neighbouring_mode_ctrls.enabled = rtc_tune ? 0 : 1;

        // near_count_ctrls
        cand_reduction_ctrls->near_count_ctrls.enabled         = 1;
        cand_reduction_ctrls->near_count_ctrls.near_count      = 1;
        cand_reduction_ctrls->near_count_ctrls.near_near_count = 1;

        // lpd1_mvp_best_me_list (LPD1 only signal)
        cand_reduction_ctrls->lpd1_mvp_best_me_list = 0;

        // cand_elimination_ctrls
        cand_reduction_ctrls->cand_elimination_ctrls.enabled         = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.dc_only         = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_me   = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_pme  = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_warp = 0;
        cand_reduction_ctrls->cand_elimination_ctrls.th_multiplier   = 1;

        // reduce_unipred_candidates
        cand_reduction_ctrls->reduce_unipred_candidates = 1;

        break;

    case 6:

        // angular reduction at mds0
        cand_reduction_ctrls->mds0_reduce_intra = 1;

        // merge_inter_classes
        cand_reduction_ctrls->merge_inter_classes = 1;

        // redundant_cand_level
        cand_reduction_ctrls->redundant_cand_ctrls.score_th = 8;
        cand_reduction_ctrls->redundant_cand_ctrls.mag_th   = 64;

        // use_neighbouring_mode
        cand_reduction_ctrls->use_neighbouring_mode_ctrls.enabled = rtc_tune ? 0 : 1;

        // near_count_ctrls
        cand_reduction_ctrls->near_count_ctrls.enabled         = 1;
        cand_reduction_ctrls->near_count_ctrls.near_count      = 1;
        cand_reduction_ctrls->near_count_ctrls.near_near_count = 1;

        // lpd1_mvp_best_me_list (LPD1 only signal)
        cand_reduction_ctrls->lpd1_mvp_best_me_list = 1;

        // cand_elimination_ctrls
        cand_reduction_ctrls->cand_elimination_ctrls.enabled         = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.dc_only         = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_me   = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_pme  = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_warp = 0;
        cand_reduction_ctrls->cand_elimination_ctrls.th_multiplier   = 1;

        // reduce_unipred_candidates
        cand_reduction_ctrls->reduce_unipred_candidates = (!pcs->ppcs->is_ref ||
                                                           ((l0_was_skip && l1_was_skip && ref_skip_perc > 35) &&
                                                            me_8x8_cost_variance < (500 * picture_qp) &&
                                                            me_64x64_distortion < (500 * picture_qp)))
            ? 3
            : 1;

        break;

    case 7:

        // angular reduction at mds0
        cand_reduction_ctrls->mds0_reduce_intra = 1;

        // merge_inter_classes
        cand_reduction_ctrls->merge_inter_classes = 1;

        // redundant_cand_level
        cand_reduction_ctrls->redundant_cand_ctrls.score_th = 8;
        cand_reduction_ctrls->redundant_cand_ctrls.mag_th   = 64;

        // use_neighbouring_mode
        cand_reduction_ctrls->use_neighbouring_mode_ctrls.enabled = rtc_tune ? 0 : 1;

        // near_count_ctrls
        cand_reduction_ctrls->near_count_ctrls.enabled         = 1;
        cand_reduction_ctrls->near_count_ctrls.near_count      = 0;
        cand_reduction_ctrls->near_count_ctrls.near_near_count = 1;

        // lpd1_mvp_best_me_list (LPD1 only signal)
        cand_reduction_ctrls->lpd1_mvp_best_me_list = 1;

        // cand_elimination_ctrls
        cand_reduction_ctrls->cand_elimination_ctrls.enabled         = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.dc_only         = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_me   = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_pme  = 1;
        cand_reduction_ctrls->cand_elimination_ctrls.inject_new_warp = 0;
        cand_reduction_ctrls->cand_elimination_ctrls.th_multiplier   = 1;

        // reduce_unipred_candidates
        cand_reduction_ctrls->reduce_unipred_candidates = (!pcs->ppcs->is_ref ||
                                                           ((l0_was_skip && l1_was_skip && ref_skip_perc > 35) &&
                                                            me_8x8_cost_variance < (500 * picture_qp) &&
                                                            me_64x64_distortion < (500 * picture_qp)))
            ? 3
            : 1;

        break;

    default: assert(0); break;
    }

    // lpd1_mvp_best_me_list can only use this feature when a single unipred ME candidate is selected,
    if (!(pcs->ppcs->ref_list0_count_try == 1 && pcs->ppcs->ref_list1_count_try == 1 &&
          pcs->ppcs->use_best_me_unipred_cand_only))
        cand_reduction_ctrls->lpd1_mvp_best_me_list = 0;
}
uint8_t svt_aom_set_chroma_controls(ModeDecisionContext *ctx, uint8_t uv_level) {
    UvCtrls *uv_ctrls = ctx ? &ctx->uv_ctrls : NULL;
    uint8_t  uv_mode  = 0;

    switch (uv_level) {
    case 0:
        uv_mode = CHROMA_MODE_2;
        if (uv_ctrls) {
            uv_ctrls->enabled = 0;
        }
        break;
    case 1:
        uv_mode = CHROMA_MODE_0;
        if (uv_ctrls) {
            uv_ctrls->enabled                = 1;
            uv_ctrls->ind_uv_last_mds        = 0;
            uv_ctrls->inter_vs_intra_cost_th = 0;
            uv_ctrls->skip_ind_uv_if_only_dc = 0;
            uv_ctrls->uv_nic_scaling_num     = 16;
        }
        break;
    case 2:
        uv_mode = CHROMA_MODE_0;
        if (uv_ctrls) {
            uv_ctrls->enabled                = 1;
            uv_ctrls->ind_uv_last_mds        = 1;
            uv_ctrls->inter_vs_intra_cost_th = 130;
            uv_ctrls->skip_ind_uv_if_only_dc = 0;
            uv_ctrls->uv_nic_scaling_num     = 8;
        }
        break;
    case 3:
        uv_mode = CHROMA_MODE_0;
        if (uv_ctrls) {
            uv_ctrls->enabled                = 1;
            uv_ctrls->ind_uv_last_mds        = 1;
            uv_ctrls->inter_vs_intra_cost_th = 100;
            uv_ctrls->skip_ind_uv_if_only_dc = 0;
            uv_ctrls->uv_nic_scaling_num     = 1;
        }
        break;
    case 4:
        uv_mode = CHROMA_MODE_0;
        if (uv_ctrls) {
            uv_ctrls->enabled                = 1;
            uv_ctrls->ind_uv_last_mds        = 2;
            uv_ctrls->inter_vs_intra_cost_th = 100;
            uv_ctrls->skip_ind_uv_if_only_dc = 1;
            uv_ctrls->uv_nic_scaling_num     = 1;
        }
        break;
    case 5:
        uv_mode = CHROMA_MODE_1;
        if (uv_ctrls) {
            uv_ctrls->enabled = 1;
        }
        break;
    default: assert(0); break;
    }

    if (ctx) {
        uv_ctrls->uv_mode = uv_mode;
    }

    return uv_mode;
}

void svt_aom_set_wm_controls(ModeDecisionContext *ctx, uint8_t wm_level) {
    WmCtrls *wm_ctrls = &ctx->wm_ctrls;

    switch (wm_level) {
    case 0: wm_ctrls->enabled = 0; break;
    case 1:
        wm_ctrls->enabled                 = 1;
        wm_ctrls->use_wm_for_mvp          = 1;
        wm_ctrls->refinement_iterations   = 8;
        wm_ctrls->refine_level            = 0;
        wm_ctrls->min_neighbour_perc      = 0;
        wm_ctrls->corner_perc_bias        = 0;
        wm_ctrls->lower_band_th           = 0;
        wm_ctrls->upper_band_th           = (uint16_t)~0;
        wm_ctrls->shut_approx_if_not_mds0 = 0;
        break;
    case 2:
        wm_ctrls->enabled                 = 1;
        wm_ctrls->use_wm_for_mvp          = 1;
        wm_ctrls->refinement_iterations   = 8;
        wm_ctrls->refine_level            = 1;
        wm_ctrls->min_neighbour_perc      = 0;
        wm_ctrls->corner_perc_bias        = 0;
        wm_ctrls->lower_band_th           = 0;
        wm_ctrls->upper_band_th           = (uint16_t)~0;
        wm_ctrls->shut_approx_if_not_mds0 = 1;
        break;
    case 3:
        wm_ctrls->enabled                 = 1;
        wm_ctrls->use_wm_for_mvp          = 1;
        wm_ctrls->refinement_iterations   = 8;
        wm_ctrls->refine_level            = 1;
        wm_ctrls->min_neighbour_perc      = 0;
        wm_ctrls->corner_perc_bias        = 0;
        wm_ctrls->lower_band_th           = 1 << 10;
        wm_ctrls->upper_band_th           = (uint16_t)~0;
        wm_ctrls->shut_approx_if_not_mds0 = 1;
        break;
    case 4:
        wm_ctrls->enabled                 = 1;
        wm_ctrls->use_wm_for_mvp          = 0;
        wm_ctrls->refinement_iterations   = 0;
        wm_ctrls->refine_level            = 0;
        wm_ctrls->min_neighbour_perc      = 0;
        wm_ctrls->corner_perc_bias        = 0;
        wm_ctrls->lower_band_th           = 0;
        wm_ctrls->upper_band_th           = (uint16_t)~0;
        wm_ctrls->shut_approx_if_not_mds0 = 0;
        break;
    default: assert(0); break;
    }
}
// Get the nic_level used for each preset (to be passed to setting function: svt_aom_set_nic_controls())
// hierarchical_levels should be the sequence-level hierarchical structure (found in scs->static_config.hierarchical_levels
uint8_t svt_aom_get_nic_level(EncMode enc_mode, uint8_t is_base, uint8_t hierarchical_levels, bool rtc_tune) {
    uint8_t nic_level;
    if (enc_mode <= ENC_MR)
        nic_level = 1;
    else if (enc_mode <= ENC_M1)
        nic_level = is_base ? 7 : 8;
    else if (enc_mode <= ENC_M2) {
        if (hierarchical_levels <= 3)
            nic_level = 10;
        else
            nic_level = 11;
    } else if (enc_mode <= ENC_M3) {
        if (hierarchical_levels <= 3)
            nic_level = 10;
        else
            nic_level = 12;
    } else if (enc_mode <= ENC_M6)
        nic_level = 13;
    else if (enc_mode <= ENC_M8)
        nic_level = 14;
    else if (enc_mode <= ENC_M10)
        nic_level = 15;
    else if (enc_mode <= ENC_M11)
        if (rtc_tune)
            nic_level = 15;
        else
            nic_level = hierarchical_levels <= 2 ? 16 : 15;
    else
        nic_level = 16;

    return nic_level;
}
/*
* Set the NIC scaling and pruning controls.
*
* This function is used in MD to set the NIC controls and is also used at memory allocation
* to allocate the candidate buffers.  Therefore, the function returns the nic_scaling_level
* (index into MD_STAGE_NICS_SCAL_NUM array).
*
* When called at memory allocation, there is no context (it is passed as NULL) so the signals
* are not set.
*/
uint8_t svt_aom_set_nic_controls(ModeDecisionContext *ctx, uint8_t nic_level) {
    NicPruningCtrls *nic_pruning_ctrls         = ctx ? &ctx->nic_ctrls.pruning_ctrls : NULL;
    uint8_t          nic_scaling_level_lte_8x8 = 0;
    uint8_t          nic_scaling_level_gt_8x8  = 0;
    uint8_t          md_staging_mode           = MD_STAGING_MODE_0;
    switch (nic_level) {
    case 0: // MAX NIC scaling; no pruning
        // NIC scaling level
        nic_scaling_level_lte_8x8 = 0;
        nic_scaling_level_gt_8x8  = 0;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = (uint64_t)~0;
            nic_pruning_ctrls->mds2_class_th = (uint64_t)~0;
            nic_pruning_ctrls->mds3_class_th = (uint64_t)~0;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th_intra  = (uint64_t)~0;
            nic_pruning_ctrls->mds1_cand_base_th_inter  = (uint64_t)~0;
            nic_pruning_ctrls->mds1_cand_th_rank_factor = 0;
            nic_pruning_ctrls->mds2_cand_base_th        = (uint64_t)~0;
            nic_pruning_ctrls->mds2_cand_th_rank_factor = 0;
            nic_pruning_ctrls->mds2_relative_dev_th     = 0;
            nic_pruning_ctrls->mds3_cand_base_th        = (uint64_t)~0;
        }
        md_staging_mode = MD_STAGING_MODE_1;
        break;

    case 1:
        // NIC scaling level
        nic_scaling_level_lte_8x8 = 0;
        nic_scaling_level_gt_8x8  = 0;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = (uint64_t)~0;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 4;

            nic_pruning_ctrls->mds3_class_th = 25;
            nic_pruning_ctrls->mds3_band_cnt = 4;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th_intra  = (uint64_t)~0;
            nic_pruning_ctrls->mds1_cand_base_th_inter  = (uint64_t)~0;
            nic_pruning_ctrls->mds1_cand_th_rank_factor = 0;
            nic_pruning_ctrls->mds2_cand_base_th        = 50;
            nic_pruning_ctrls->mds2_cand_th_rank_factor = 0;
            nic_pruning_ctrls->mds2_relative_dev_th     = 0;
            nic_pruning_ctrls->mds3_cand_base_th        = 50;
        }
        md_staging_mode = MD_STAGING_MODE_2;
        break;

    case 2:
        // NIC scaling level
        nic_scaling_level_lte_8x8 = 1;
        nic_scaling_level_gt_8x8  = 1;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = (uint64_t)~0;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 4;

            nic_pruning_ctrls->mds3_class_th = 25;
            nic_pruning_ctrls->mds3_band_cnt = 4;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th_intra  = (uint64_t)~0;
            nic_pruning_ctrls->mds1_cand_base_th_inter  = (uint64_t)~0;
            nic_pruning_ctrls->mds1_cand_th_rank_factor = 0;
            nic_pruning_ctrls->mds2_cand_base_th        = 50;
            nic_pruning_ctrls->mds2_cand_th_rank_factor = 0;
            nic_pruning_ctrls->mds2_relative_dev_th     = 0;
            nic_pruning_ctrls->mds3_cand_base_th        = 50;
        }
        md_staging_mode = MD_STAGING_MODE_2;
        break;

    case 3:
        // NIC scaling level
        nic_scaling_level_lte_8x8 = 1;
        nic_scaling_level_gt_8x8  = 1;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = (uint64_t)~0;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 4;

            nic_pruning_ctrls->mds3_class_th = 25;
            nic_pruning_ctrls->mds3_band_cnt = 8;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th_intra  = 1200;
            nic_pruning_ctrls->mds1_cand_base_th_inter  = 500;
            nic_pruning_ctrls->mds1_cand_th_rank_factor = 0;
            nic_pruning_ctrls->mds2_cand_base_th        = 30;
            nic_pruning_ctrls->mds2_cand_th_rank_factor = 0;
            nic_pruning_ctrls->mds2_relative_dev_th     = 0;
            nic_pruning_ctrls->mds3_cand_base_th        = 30;
        }
        md_staging_mode = MD_STAGING_MODE_2;
        break;

    case 4:
        // NIC scaling level
        nic_scaling_level_lte_8x8 = 2;
        nic_scaling_level_gt_8x8  = 2;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = (uint64_t)~0;
            nic_pruning_ctrls->mds1_band_cnt = 2;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 4;

            nic_pruning_ctrls->mds3_class_th = 25;
            nic_pruning_ctrls->mds3_band_cnt = 12;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th_intra  = 1200;
            nic_pruning_ctrls->mds1_cand_base_th_inter  = 300;
            nic_pruning_ctrls->mds1_cand_th_rank_factor = 0;
            nic_pruning_ctrls->mds2_cand_base_th        = 20;
            nic_pruning_ctrls->mds2_cand_th_rank_factor = 0;
            nic_pruning_ctrls->mds2_relative_dev_th     = 0;
            nic_pruning_ctrls->mds3_cand_base_th        = 20;
        }
        md_staging_mode = MD_STAGING_MODE_2;
        break;

    case 5:
        // NIC scaling level
        nic_scaling_level_lte_8x8 = 2;
        nic_scaling_level_gt_8x8  = 2;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 500;
            nic_pruning_ctrls->mds1_band_cnt = 3;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 8;

            nic_pruning_ctrls->mds3_class_th = 20;
            nic_pruning_ctrls->mds3_band_cnt = 12;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th_intra  = 1200;
            nic_pruning_ctrls->mds1_cand_base_th_inter  = 300;
            nic_pruning_ctrls->mds1_cand_th_rank_factor = 0;
            nic_pruning_ctrls->mds2_cand_base_th        = 20;
            nic_pruning_ctrls->mds2_cand_th_rank_factor = 0;
            nic_pruning_ctrls->mds2_relative_dev_th     = 0;
            nic_pruning_ctrls->mds3_cand_base_th        = 15;
        }
        md_staging_mode = MD_STAGING_MODE_2;
        break;

    case 6:
        // NIC scaling level
        nic_scaling_level_lte_8x8 = 3;
        nic_scaling_level_gt_8x8  = 3;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 500;
            nic_pruning_ctrls->mds1_band_cnt = 3;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 8;

            nic_pruning_ctrls->mds3_class_th = 20;
            nic_pruning_ctrls->mds3_band_cnt = 12;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th_intra  = 1200;
            nic_pruning_ctrls->mds1_cand_base_th_inter  = 300;
            nic_pruning_ctrls->mds1_cand_th_rank_factor = 0;
            nic_pruning_ctrls->mds2_cand_base_th        = 20;
            nic_pruning_ctrls->mds2_cand_th_rank_factor = 0;
            nic_pruning_ctrls->mds2_relative_dev_th     = 0;
            nic_pruning_ctrls->mds3_cand_base_th        = 15;
        }
        md_staging_mode = MD_STAGING_MODE_2;
        break;

    case 7:
        // NIC scaling level
        nic_scaling_level_lte_8x8 = 4;
        nic_scaling_level_gt_8x8  = 4;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 500;
            nic_pruning_ctrls->mds1_band_cnt = 3;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 8;

            nic_pruning_ctrls->mds3_class_th = 20;
            nic_pruning_ctrls->mds3_band_cnt = 12;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th_intra  = 1200;
            nic_pruning_ctrls->mds1_cand_base_th_inter  = 300;
            nic_pruning_ctrls->mds1_cand_th_rank_factor = 0;
            nic_pruning_ctrls->mds2_cand_base_th        = 20;
            nic_pruning_ctrls->mds2_cand_th_rank_factor = 0;
            nic_pruning_ctrls->mds2_relative_dev_th     = 0;
            nic_pruning_ctrls->mds3_cand_base_th        = 15;
        }
        md_staging_mode = MD_STAGING_MODE_2;
        break;
    case 8:
        // NIC scaling level
        nic_scaling_level_lte_8x8 = 6;
        nic_scaling_level_gt_8x8  = 6;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 500;
            nic_pruning_ctrls->mds1_band_cnt = 3;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 8;

            nic_pruning_ctrls->mds3_class_th = 20;
            nic_pruning_ctrls->mds3_band_cnt = 12;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th_intra  = 1200;
            nic_pruning_ctrls->mds1_cand_base_th_inter  = 300;
            nic_pruning_ctrls->mds1_cand_th_rank_factor = 0;
            nic_pruning_ctrls->mds2_cand_base_th        = 20;
            nic_pruning_ctrls->mds2_cand_th_rank_factor = 0;
            nic_pruning_ctrls->mds2_relative_dev_th     = 0;
            nic_pruning_ctrls->mds3_cand_base_th        = 15;
        }
        md_staging_mode = MD_STAGING_MODE_2;
        break;
    case 9:
        // NIC scaling level
        nic_scaling_level_lte_8x8 = 6;
        nic_scaling_level_gt_8x8  = 6;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 300;
            nic_pruning_ctrls->mds1_band_cnt = 4;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 10;

            nic_pruning_ctrls->mds3_class_th = 15;
            nic_pruning_ctrls->mds3_band_cnt = 16;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th_intra  = 1200;
            nic_pruning_ctrls->mds1_cand_base_th_inter  = 300;
            nic_pruning_ctrls->mds1_cand_th_rank_factor = 0;
            nic_pruning_ctrls->mds2_cand_base_th        = 20;
            nic_pruning_ctrls->mds2_cand_th_rank_factor = 0;
            nic_pruning_ctrls->mds2_relative_dev_th     = 0;
            nic_pruning_ctrls->mds3_cand_base_th        = 15;
        }
        md_staging_mode = MD_STAGING_MODE_2;
        break;
    case 10:
        // NIC scaling level
        nic_scaling_level_lte_8x8 = 8;
        nic_scaling_level_gt_8x8  = 8;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 300;
            nic_pruning_ctrls->mds1_band_cnt = 4;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 10;

            nic_pruning_ctrls->mds3_class_th = 15;
            nic_pruning_ctrls->mds3_band_cnt = 16;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th_intra  = 1200;
            nic_pruning_ctrls->mds1_cand_base_th_inter  = 300;
            nic_pruning_ctrls->mds1_cand_th_rank_factor = 0;
            nic_pruning_ctrls->mds2_cand_base_th        = 20;
            nic_pruning_ctrls->mds2_cand_th_rank_factor = 0;
            nic_pruning_ctrls->mds2_relative_dev_th     = 0;
            nic_pruning_ctrls->mds3_cand_base_th        = 15;
        }
        md_staging_mode = MD_STAGING_MODE_2;
        break;
    case 11:
        // NIC scaling level
        nic_scaling_level_lte_8x8 = 9;
        nic_scaling_level_gt_8x8  = 10;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 300;
            nic_pruning_ctrls->mds1_band_cnt = 4;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 10;

            nic_pruning_ctrls->mds3_class_th = 15;
            nic_pruning_ctrls->mds3_band_cnt = 16;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th_intra  = 1200;
            nic_pruning_ctrls->mds1_cand_base_th_inter  = 300;
            nic_pruning_ctrls->mds1_cand_th_rank_factor = 0;

            nic_pruning_ctrls->mds2_cand_base_th        = 20;
            nic_pruning_ctrls->mds2_cand_th_rank_factor = 1;
            nic_pruning_ctrls->mds2_relative_dev_th     = 5;

            nic_pruning_ctrls->mds3_cand_base_th = 15;
        }
        md_staging_mode = MD_STAGING_MODE_1;
        break;
    case 12:
        // NIC scaling level
        nic_scaling_level_lte_8x8 = 9;
        nic_scaling_level_gt_8x8  = 10;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 300;
            nic_pruning_ctrls->mds1_band_cnt = 4;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 10;

            nic_pruning_ctrls->mds3_class_th = 15;
            nic_pruning_ctrls->mds3_band_cnt = 16;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th_intra  = 1200;
            nic_pruning_ctrls->mds1_cand_base_th_inter  = 300;
            nic_pruning_ctrls->mds1_cand_th_rank_factor = 3;
            nic_pruning_ctrls->mds2_cand_base_th        = 15;
            nic_pruning_ctrls->mds2_cand_th_rank_factor = 1;
            nic_pruning_ctrls->mds2_relative_dev_th     = 5;
            nic_pruning_ctrls->mds3_cand_base_th        = 15;
        }
        md_staging_mode = MD_STAGING_MODE_1;
        break;
    case 13:
        // NIC scaling level
        nic_scaling_level_lte_8x8 = 9;
        nic_scaling_level_gt_8x8  = 12;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 200;
            nic_pruning_ctrls->mds1_band_cnt = 16;

            nic_pruning_ctrls->mds2_class_th = 25;
            nic_pruning_ctrls->mds2_band_cnt = 10;

            nic_pruning_ctrls->mds3_class_th = 15;
            nic_pruning_ctrls->mds3_band_cnt = 16;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th_intra  = 1200;
            nic_pruning_ctrls->mds1_cand_base_th_inter  = 300;
            nic_pruning_ctrls->mds1_cand_th_rank_factor = 3;
            nic_pruning_ctrls->mds2_cand_base_th        = 15;
            nic_pruning_ctrls->mds2_cand_th_rank_factor = 1;
            nic_pruning_ctrls->mds2_relative_dev_th     = 5;
            nic_pruning_ctrls->mds3_cand_base_th        = 15;
        }
        md_staging_mode = MD_STAGING_MODE_1;
        break;
    case 14:
        // NIC scaling level
        nic_scaling_level_lte_8x8 = 13;
        nic_scaling_level_gt_8x8  = 14;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 200;
            nic_pruning_ctrls->mds1_band_cnt = 16;

            nic_pruning_ctrls->mds2_class_th = 10;
            nic_pruning_ctrls->mds2_band_cnt = 2;

            nic_pruning_ctrls->mds3_class_th = 10;
            nic_pruning_ctrls->mds3_band_cnt = 16;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th_intra  = 300;
            nic_pruning_ctrls->mds1_cand_base_th_inter  = 300;
            nic_pruning_ctrls->mds1_cand_th_rank_factor = 3;
            nic_pruning_ctrls->mds2_cand_base_th        = 5;
            nic_pruning_ctrls->mds2_cand_th_rank_factor = 1;
            nic_pruning_ctrls->mds2_relative_dev_th     = 5;
            nic_pruning_ctrls->mds3_cand_base_th        = 5;
        }
        md_staging_mode = MD_STAGING_MODE_1;
        break;
    case 15:
        // NIC scaling level
        nic_scaling_level_lte_8x8 = 14;
        nic_scaling_level_gt_8x8  = 14;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 200;
            nic_pruning_ctrls->mds1_band_cnt = 16;

            nic_pruning_ctrls->mds2_class_th = 10;
            nic_pruning_ctrls->mds2_band_cnt = 2;

            nic_pruning_ctrls->mds3_class_th = 10;
            nic_pruning_ctrls->mds3_band_cnt = 16;

            nic_pruning_ctrls->enable_skipping_mds1 = 0;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th_intra  = 300;
            nic_pruning_ctrls->mds1_cand_base_th_inter  = 300;
            nic_pruning_ctrls->mds1_cand_th_rank_factor = 3;
            nic_pruning_ctrls->mds2_cand_base_th        = 5;
            nic_pruning_ctrls->mds2_cand_th_rank_factor = 1;
            nic_pruning_ctrls->mds2_relative_dev_th     = 5;
            nic_pruning_ctrls->mds3_cand_base_th        = 5;
        }
        md_staging_mode = MD_STAGING_MODE_1;
        break;

    case 16:
        // NIC scaling level
        nic_scaling_level_lte_8x8 = 15;
        nic_scaling_level_gt_8x8  = 15;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 75;
            nic_pruning_ctrls->mds1_band_cnt = 16;

            nic_pruning_ctrls->mds2_class_th = 0;
            nic_pruning_ctrls->mds2_band_cnt = 2;

            nic_pruning_ctrls->mds3_class_th = 0;
            nic_pruning_ctrls->mds3_band_cnt = 2;

            nic_pruning_ctrls->enable_skipping_mds1 = 1;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th_intra  = 1;
            nic_pruning_ctrls->mds1_cand_base_th_inter  = 1;
            nic_pruning_ctrls->mds1_cand_th_rank_factor = 3;
            nic_pruning_ctrls->mds2_cand_base_th        = 1;
            nic_pruning_ctrls->mds2_cand_th_rank_factor = 1;
            nic_pruning_ctrls->mds2_relative_dev_th     = 5;
            nic_pruning_ctrls->mds3_cand_base_th        = 1;
        }
        md_staging_mode = MD_STAGING_MODE_1;
        break;

    case 17:
        // NIC scaling level
        nic_scaling_level_lte_8x8 = 15;
        nic_scaling_level_gt_8x8  = 15;

        if (nic_pruning_ctrls) {
            // Class pruning settings
            nic_pruning_ctrls->mds1_class_th = 75;
            nic_pruning_ctrls->mds1_band_cnt = 16;

            nic_pruning_ctrls->mds2_class_th = 0;
            nic_pruning_ctrls->mds2_band_cnt = 2;

            nic_pruning_ctrls->mds3_class_th = 0;
            nic_pruning_ctrls->mds3_band_cnt = 2;

            nic_pruning_ctrls->enable_skipping_mds1 = 1;
            nic_pruning_ctrls->force_1_cand_th      = 0;

            // Cand pruning settings
            nic_pruning_ctrls->mds1_cand_base_th_intra  = 1;
            nic_pruning_ctrls->mds1_cand_base_th_inter  = 1;
            nic_pruning_ctrls->mds1_cand_th_rank_factor = 3;
            nic_pruning_ctrls->mds2_cand_base_th        = 1;
            nic_pruning_ctrls->mds2_cand_th_rank_factor = 1;
            nic_pruning_ctrls->mds2_relative_dev_th     = 5;
            nic_pruning_ctrls->mds3_cand_base_th        = 1;
        }
        md_staging_mode = MD_STAGING_MODE_0;
        break;
    default: assert(0); break;
    }

    if (ctx) {
        uint8_t nic_scaling_level = ctx->blk_geom->sq_size <= 8 ? nic_scaling_level_lte_8x8 : nic_scaling_level_gt_8x8;
        NicScalingCtrls *nic_scaling_ctrls    = &ctx->nic_ctrls.scaling_ctrls;
        nic_scaling_ctrls->stage1_scaling_num = MD_STAGE_NICS_SCAL_NUM[nic_scaling_level][MD_STAGE_1];
        nic_scaling_ctrls->stage2_scaling_num = MD_STAGE_NICS_SCAL_NUM[nic_scaling_level][MD_STAGE_2];
        nic_scaling_ctrls->stage3_scaling_num = MD_STAGE_NICS_SCAL_NUM[nic_scaling_level][MD_STAGE_3];
        ctx->nic_ctrls.md_staging_mode        = md_staging_mode;
    }
    // return the lowest level of NIC scaling that can be used for memory allocation
    uint8_t nic_scaling_level = MIN(nic_scaling_level_lte_8x8, nic_scaling_level_gt_8x8);
    return nic_scaling_level;
}
/*
* This function is used in MD to set the NSQ controls.
*/
void svt_aom_set_nsq_ctrls(ModeDecisionContext *ctx, uint8_t nsq_level, uint8_t *allow_HVA_HVB, uint8_t *allow_HV4,
                           uint8_t *min_nsq_bsize) {
    NsqCtrls  nsq_ctrls_struct;
    NsqCtrls *nsq_ctrls = &nsq_ctrls_struct;
    switch (nsq_level) {
    case 0:
        nsq_ctrls->enabled               = 0;
        nsq_ctrls->psq_cplx_lvl          = 0;
        nsq_ctrls->allow_HV4             = 0;
        nsq_ctrls->allow_HVA_HVB         = 0;
        nsq_ctrls->min_nsq_block_size    = 0;
        nsq_ctrls->psq_txs_lvl           = 0;
        nsq_ctrls->sub_depth_block_lvl   = 0;
        nsq_ctrls->psq_pred_lvl          = 0;
        nsq_ctrls->component_multiple_th = 0;
        break;

    case 1: // Original MRS level
        nsq_ctrls->enabled = 1;

        nsq_ctrls->min_nsq_block_size = 0;
        nsq_ctrls->allow_HV4          = 1;
        nsq_ctrls->allow_HVA_HVB      = 1;

        nsq_ctrls->sq_weight                    = (uint32_t)~0;
        nsq_ctrls->psq_cplx_lvl                 = 0;
        nsq_ctrls->max_part0_to_part1_dev       = 0;
        nsq_ctrls->nsq_split_cost_th            = 0;
        nsq_ctrls->lower_depth_split_cost_th    = 0;
        nsq_ctrls->H_vs_V_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_modulation = 0;
        nsq_ctrls->psq_txs_lvl                  = 0;
        nsq_ctrls->sub_depth_block_lvl          = 0;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 0;
        nsq_ctrls->hv_weight                    = 110;
        break;

    case 2:
        nsq_ctrls->enabled = 1;

        nsq_ctrls->min_nsq_block_size = 0;
        nsq_ctrls->allow_HV4          = 1;
        nsq_ctrls->allow_HVA_HVB      = 1;

        nsq_ctrls->sq_weight                    = 105;
        nsq_ctrls->psq_cplx_lvl                 = 0;
        nsq_ctrls->max_part0_to_part1_dev       = 0;
        nsq_ctrls->nsq_split_cost_th            = 150;
        nsq_ctrls->lower_depth_split_cost_th    = 5;
        nsq_ctrls->H_vs_V_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_modulation = 0;
        nsq_ctrls->psq_txs_lvl                  = 0;
        nsq_ctrls->sub_depth_block_lvl          = 0;

        nsq_ctrls->psq_pred_lvl          = 0;
        nsq_ctrls->component_multiple_th = 0;
        nsq_ctrls->hv_weight             = 110;
        break;

    case 3:
        nsq_ctrls->enabled = 1;

        nsq_ctrls->min_nsq_block_size = 0;
        nsq_ctrls->allow_HV4          = 1;
        nsq_ctrls->allow_HVA_HVB      = 1;

        nsq_ctrls->sq_weight                    = 105;
        nsq_ctrls->psq_cplx_lvl                 = 1;
        nsq_ctrls->max_part0_to_part1_dev       = 0;
        nsq_ctrls->nsq_split_cost_th            = 150;
        nsq_ctrls->lower_depth_split_cost_th    = 5;
        nsq_ctrls->H_vs_V_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_modulation = 0;
        nsq_ctrls->psq_txs_lvl                  = 0;
        nsq_ctrls->sub_depth_block_lvl          = 0;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 0;
        nsq_ctrls->hv_weight                    = 110;
        break;

    case 4:
        nsq_ctrls->enabled = 1;

        nsq_ctrls->min_nsq_block_size = 0;
        nsq_ctrls->allow_HV4          = 1;
        nsq_ctrls->allow_HVA_HVB      = 1;

        nsq_ctrls->sq_weight                    = 105;
        nsq_ctrls->psq_cplx_lvl                 = 2;
        nsq_ctrls->max_part0_to_part1_dev       = 0;
        nsq_ctrls->nsq_split_cost_th            = 150;
        nsq_ctrls->lower_depth_split_cost_th    = 5;
        nsq_ctrls->H_vs_V_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_modulation = 0;
        nsq_ctrls->psq_txs_lvl                  = 0;
        nsq_ctrls->sub_depth_block_lvl          = 0;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 0;
        nsq_ctrls->hv_weight                    = 110;
        break;

    case 5:
        nsq_ctrls->enabled = 1;

        nsq_ctrls->min_nsq_block_size = 0;
        nsq_ctrls->allow_HV4          = 1;
        nsq_ctrls->allow_HVA_HVB      = 0;

        nsq_ctrls->sq_weight                    = 100;
        nsq_ctrls->psq_cplx_lvl                 = 2;
        nsq_ctrls->max_part0_to_part1_dev       = 0;
        nsq_ctrls->nsq_split_cost_th            = 100;
        nsq_ctrls->lower_depth_split_cost_th    = 5;
        nsq_ctrls->H_vs_V_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_modulation = 0;
        nsq_ctrls->psq_txs_lvl                  = 0;
        nsq_ctrls->sub_depth_block_lvl          = 1;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 80;
        nsq_ctrls->hv_weight                    = 110;
        break;

    case 6:
        nsq_ctrls->enabled = 1;

        nsq_ctrls->min_nsq_block_size = 0;
        nsq_ctrls->allow_HV4          = 1;
        nsq_ctrls->allow_HVA_HVB      = 0;

        nsq_ctrls->sq_weight                    = 100;
        nsq_ctrls->psq_cplx_lvl                 = 3;
        nsq_ctrls->max_part0_to_part1_dev       = 0;
        nsq_ctrls->nsq_split_cost_th            = 100;
        nsq_ctrls->lower_depth_split_cost_th    = 5;
        nsq_ctrls->H_vs_V_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_modulation = 0;
        nsq_ctrls->psq_txs_lvl                  = 0;
        nsq_ctrls->sub_depth_block_lvl          = 1;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 80;
        nsq_ctrls->hv_weight                    = 110;
        break;

    case 7:
        nsq_ctrls->enabled = 1;

        nsq_ctrls->min_nsq_block_size = 0;
        nsq_ctrls->allow_HV4          = 1;
        nsq_ctrls->allow_HVA_HVB      = 0;

        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 2;
        nsq_ctrls->max_part0_to_part1_dev       = 0;
        nsq_ctrls->nsq_split_cost_th            = 80;
        nsq_ctrls->lower_depth_split_cost_th    = 5;
        nsq_ctrls->H_vs_V_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_modulation = 0;
        nsq_ctrls->psq_txs_lvl                  = 0;
        nsq_ctrls->sub_depth_block_lvl          = 1;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 80;
        nsq_ctrls->hv_weight                    = 100;
        break;

    case 8:
        nsq_ctrls->enabled = 1;

        nsq_ctrls->min_nsq_block_size = 0;
        nsq_ctrls->allow_HV4          = 1;
        nsq_ctrls->allow_HVA_HVB      = 0;

        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 3;
        nsq_ctrls->max_part0_to_part1_dev       = 0;
        nsq_ctrls->nsq_split_cost_th            = 80;
        nsq_ctrls->lower_depth_split_cost_th    = 5;
        nsq_ctrls->H_vs_V_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_modulation = 0;
        nsq_ctrls->psq_txs_lvl                  = 0;
        nsq_ctrls->sub_depth_block_lvl          = 1;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 80;
        nsq_ctrls->hv_weight                    = 100;
        break;
    case 9:
        nsq_ctrls->enabled = 1;

        nsq_ctrls->min_nsq_block_size = 0;
        nsq_ctrls->allow_HV4          = 1;
        nsq_ctrls->allow_HVA_HVB      = 0;

        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 2;
        nsq_ctrls->max_part0_to_part1_dev       = 0;
        nsq_ctrls->nsq_split_cost_th            = 80;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_modulation = 0;
        nsq_ctrls->psq_txs_lvl                  = 0;
        nsq_ctrls->sub_depth_block_lvl          = 1;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 60;
        nsq_ctrls->hv_weight                    = 100;
        break;

    case 10:
        nsq_ctrls->enabled = 1;

        nsq_ctrls->min_nsq_block_size = 0;
        nsq_ctrls->allow_HV4          = 1;
        nsq_ctrls->allow_HVA_HVB      = 0;

        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 3;
        nsq_ctrls->max_part0_to_part1_dev       = 0;
        nsq_ctrls->nsq_split_cost_th            = 80;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_modulation = 0;
        nsq_ctrls->psq_txs_lvl                  = 0;
        nsq_ctrls->sub_depth_block_lvl          = 1;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 60;
        nsq_ctrls->hv_weight                    = 100;
        break;

    case 11:
        nsq_ctrls->enabled = 1;

        nsq_ctrls->min_nsq_block_size = 8;
        nsq_ctrls->allow_HV4          = 1;
        nsq_ctrls->allow_HVA_HVB      = 0;

        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 3;
        nsq_ctrls->max_part0_to_part1_dev       = 0;
        nsq_ctrls->nsq_split_cost_th            = 80;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_th         = 10;
        nsq_ctrls->non_HV_split_rate_modulation = 1;
        nsq_ctrls->psq_txs_lvl                  = 0;
        nsq_ctrls->sub_depth_block_lvl          = 1;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 20;
        nsq_ctrls->hv_weight                    = 100;
        break;

    case 12:
        nsq_ctrls->enabled = 1;

        nsq_ctrls->min_nsq_block_size = 8;
        nsq_ctrls->allow_HV4          = 1;
        nsq_ctrls->allow_HVA_HVB      = 0;

        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 3;
        nsq_ctrls->max_part0_to_part1_dev       = 40;
        nsq_ctrls->nsq_split_cost_th            = 80;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_th         = 10;
        nsq_ctrls->non_HV_split_rate_modulation = 1;
        nsq_ctrls->psq_txs_lvl                  = 0;
        nsq_ctrls->sub_depth_block_lvl          = 1;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 20;
        nsq_ctrls->hv_weight                    = 100;
        break;

    case 13:
        nsq_ctrls->enabled = 1;

        nsq_ctrls->min_nsq_block_size = 8;
        nsq_ctrls->allow_HV4          = 1;
        nsq_ctrls->allow_HVA_HVB      = 0;

        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 3;
        nsq_ctrls->max_part0_to_part1_dev       = 60;
        nsq_ctrls->nsq_split_cost_th            = 80;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_th         = 10;
        nsq_ctrls->non_HV_split_rate_modulation = 1;
        nsq_ctrls->psq_txs_lvl                  = 0;
        nsq_ctrls->sub_depth_block_lvl          = 1;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 20;
        nsq_ctrls->hv_weight                    = 100;
        break;

    case 14:
        nsq_ctrls->enabled = 1;

        nsq_ctrls->min_nsq_block_size = 8;
        nsq_ctrls->allow_HV4          = 1;
        nsq_ctrls->allow_HVA_HVB      = 0;

        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 3;
        nsq_ctrls->max_part0_to_part1_dev       = 0;
        nsq_ctrls->nsq_split_cost_th            = 80;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_th         = 10;
        nsq_ctrls->non_HV_split_rate_modulation = 1;
        nsq_ctrls->psq_txs_lvl                  = 0;
        nsq_ctrls->sub_depth_block_lvl          = 1;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 20;
        nsq_ctrls->hv_weight                    = 100;
        break;

    case 15:
        nsq_ctrls->enabled = 1;

        nsq_ctrls->min_nsq_block_size = 8;
        nsq_ctrls->allow_HV4          = 1;
        nsq_ctrls->allow_HVA_HVB      = 0;

        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 3;
        nsq_ctrls->max_part0_to_part1_dev       = 40;
        nsq_ctrls->nsq_split_cost_th            = 80;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_th         = 10;
        nsq_ctrls->non_HV_split_rate_modulation = 1;
        nsq_ctrls->psq_txs_lvl                  = 0;
        nsq_ctrls->sub_depth_block_lvl          = 1;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 20;
        nsq_ctrls->hv_weight                    = 100;
        break;

    case 16:
        nsq_ctrls->enabled = 1;

        nsq_ctrls->min_nsq_block_size = 8;
        nsq_ctrls->allow_HV4          = 1;
        nsq_ctrls->allow_HVA_HVB      = 0;

        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 3;
        nsq_ctrls->max_part0_to_part1_dev       = 60;
        nsq_ctrls->nsq_split_cost_th            = 80;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_th         = 10;
        nsq_ctrls->non_HV_split_rate_modulation = 1;
        nsq_ctrls->psq_txs_lvl                  = 0;
        nsq_ctrls->sub_depth_block_lvl          = 1;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 20;
        nsq_ctrls->hv_weight                    = 100;
        break;
    case 17:
        nsq_ctrls->enabled = 1;

        nsq_ctrls->min_nsq_block_size = 8;
        nsq_ctrls->allow_HV4          = 1;
        nsq_ctrls->allow_HVA_HVB      = 0;

        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 3;
        nsq_ctrls->max_part0_to_part1_dev       = 0;
        nsq_ctrls->nsq_split_cost_th            = 60;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 60;
        nsq_ctrls->non_HV_split_rate_th         = 30;
        nsq_ctrls->non_HV_split_rate_modulation = 1;
        nsq_ctrls->psq_txs_lvl                  = 0;
        nsq_ctrls->sub_depth_block_lvl          = 1;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 20;
        nsq_ctrls->hv_weight                    = 100;
        break;
    case 18:
        nsq_ctrls->enabled = 1;

        nsq_ctrls->min_nsq_block_size = 8;
        nsq_ctrls->allow_HV4          = 1;
        nsq_ctrls->allow_HVA_HVB      = 0;

        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 3;
        nsq_ctrls->max_part0_to_part1_dev       = 40;
        nsq_ctrls->nsq_split_cost_th            = 60;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 60;
        nsq_ctrls->non_HV_split_rate_th         = 30;
        nsq_ctrls->non_HV_split_rate_modulation = 1;
        nsq_ctrls->psq_txs_lvl                  = 0;
        nsq_ctrls->sub_depth_block_lvl          = 1;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 20;
        nsq_ctrls->hv_weight                    = 100;
        break;
    case 19:
        nsq_ctrls->enabled = 1;

        nsq_ctrls->min_nsq_block_size = 8;
        nsq_ctrls->allow_HV4          = 1;
        nsq_ctrls->allow_HVA_HVB      = 0;

        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 3;
        nsq_ctrls->max_part0_to_part1_dev       = 60;
        nsq_ctrls->nsq_split_cost_th            = 60;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 60;
        nsq_ctrls->non_HV_split_rate_th         = 30;
        nsq_ctrls->non_HV_split_rate_modulation = 1;
        nsq_ctrls->psq_txs_lvl                  = 0;
        nsq_ctrls->sub_depth_block_lvl          = 1;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 20;
        nsq_ctrls->hv_weight                    = 100;
        break;
    case 20:
        nsq_ctrls->enabled = 1;

        nsq_ctrls->min_nsq_block_size = 8;
        nsq_ctrls->allow_HV4          = 1;
        nsq_ctrls->allow_HVA_HVB      = 0;

        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 4;
        nsq_ctrls->max_part0_to_part1_dev       = 0;
        nsq_ctrls->nsq_split_cost_th            = 60;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 60;
        nsq_ctrls->non_HV_split_rate_th         = 30;
        nsq_ctrls->non_HV_split_rate_modulation = 1;
        nsq_ctrls->psq_txs_lvl                  = 0;
        nsq_ctrls->sub_depth_block_lvl          = 1;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 20;
        nsq_ctrls->hv_weight                    = 100;
        break;
    case 21:
        nsq_ctrls->enabled = 1;

        nsq_ctrls->min_nsq_block_size = 8;
        nsq_ctrls->allow_HV4          = 0;
        nsq_ctrls->allow_HVA_HVB      = 0;

        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 4;
        nsq_ctrls->max_part0_to_part1_dev       = 40;
        nsq_ctrls->nsq_split_cost_th            = 60;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 60;
        nsq_ctrls->non_HV_split_rate_th         = 30;
        nsq_ctrls->non_HV_split_rate_modulation = 1;
        nsq_ctrls->psq_txs_lvl                  = 0;
        nsq_ctrls->sub_depth_block_lvl          = 1;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 20;
        nsq_ctrls->hv_weight                    = 100;
        break;
    case 22:
        nsq_ctrls->enabled            = 1;
        nsq_ctrls->min_nsq_block_size = 8;
        nsq_ctrls->allow_HV4          = 0;
        nsq_ctrls->allow_HVA_HVB      = 0;

        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 5;
        nsq_ctrls->max_part0_to_part1_dev       = 60;
        nsq_ctrls->nsq_split_cost_th            = 60;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 60;
        nsq_ctrls->non_HV_split_rate_th         = 30;
        nsq_ctrls->non_HV_split_rate_modulation = 1;
        nsq_ctrls->psq_txs_lvl                  = 0;
        nsq_ctrls->sub_depth_block_lvl          = 1;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 20;
        nsq_ctrls->hv_weight                    = 100;
        break;
    case 23:
        nsq_ctrls->enabled                      = 1;
        nsq_ctrls->min_nsq_block_size           = 8;
        nsq_ctrls->allow_HV4                    = 0;
        nsq_ctrls->allow_HVA_HVB                = 0;
        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 3;
        nsq_ctrls->max_part0_to_part1_dev       = 0;
        nsq_ctrls->nsq_split_cost_th            = 60;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 60;
        nsq_ctrls->non_HV_split_rate_th         = 30;
        nsq_ctrls->non_HV_split_rate_modulation = 1;
        nsq_ctrls->psq_txs_lvl                  = 1;
        nsq_ctrls->sub_depth_block_lvl          = 2;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 10;
        nsq_ctrls->hv_weight                    = 100;
        break;

    case 24:
        nsq_ctrls->enabled                      = 1;
        nsq_ctrls->min_nsq_block_size           = 8;
        nsq_ctrls->allow_HV4                    = 0;
        nsq_ctrls->allow_HVA_HVB                = 0;
        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 3;
        nsq_ctrls->max_part0_to_part1_dev       = 40;
        nsq_ctrls->nsq_split_cost_th            = 60;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 60;
        nsq_ctrls->non_HV_split_rate_th         = 30;
        nsq_ctrls->non_HV_split_rate_modulation = 1;
        nsq_ctrls->psq_txs_lvl                  = 1;
        nsq_ctrls->sub_depth_block_lvl          = 2;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 10;
        nsq_ctrls->hv_weight                    = 100;
        break;

    case 25:
        nsq_ctrls->enabled                      = 1;
        nsq_ctrls->min_nsq_block_size           = 8;
        nsq_ctrls->allow_HV4                    = 0;
        nsq_ctrls->allow_HVA_HVB                = 0;
        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 3;
        nsq_ctrls->max_part0_to_part1_dev       = 60;
        nsq_ctrls->nsq_split_cost_th            = 60;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 60;
        nsq_ctrls->non_HV_split_rate_th         = 30;
        nsq_ctrls->non_HV_split_rate_modulation = 1;
        nsq_ctrls->psq_txs_lvl                  = 1;
        nsq_ctrls->sub_depth_block_lvl          = 2;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 10;
        nsq_ctrls->hv_weight                    = 100;
        break;

    case 26:
        nsq_ctrls->enabled                      = 1;
        nsq_ctrls->min_nsq_block_size           = 8;
        nsq_ctrls->allow_HV4                    = 0;
        nsq_ctrls->allow_HVA_HVB                = 0;
        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 5;
        nsq_ctrls->max_part0_to_part1_dev       = 0;
        nsq_ctrls->nsq_split_cost_th            = 60;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 60;
        nsq_ctrls->non_HV_split_rate_th         = 30;
        nsq_ctrls->non_HV_split_rate_modulation = 1;
        nsq_ctrls->psq_txs_lvl                  = 1;
        nsq_ctrls->sub_depth_block_lvl          = 2;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 10;
        nsq_ctrls->hv_weight                    = 100;
        break;

    case 27:
        nsq_ctrls->enabled                      = 1;
        nsq_ctrls->min_nsq_block_size           = 8;
        nsq_ctrls->allow_HV4                    = 0;
        nsq_ctrls->allow_HVA_HVB                = 0;
        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 5;
        nsq_ctrls->max_part0_to_part1_dev       = 40;
        nsq_ctrls->nsq_split_cost_th            = 60;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 60;
        nsq_ctrls->non_HV_split_rate_th         = 30;
        nsq_ctrls->non_HV_split_rate_modulation = 1;
        nsq_ctrls->psq_txs_lvl                  = 1;
        nsq_ctrls->sub_depth_block_lvl          = 2;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 10;
        nsq_ctrls->hv_weight                    = 100;
        break;

    case 28:
        nsq_ctrls->enabled                      = 1;
        nsq_ctrls->min_nsq_block_size           = 8;
        nsq_ctrls->allow_HV4                    = 0;
        nsq_ctrls->allow_HVA_HVB                = 0;
        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 5;
        nsq_ctrls->max_part0_to_part1_dev       = 60;
        nsq_ctrls->nsq_split_cost_th            = 60;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 60;
        nsq_ctrls->non_HV_split_rate_th         = 30;
        nsq_ctrls->non_HV_split_rate_modulation = 1;
        nsq_ctrls->psq_txs_lvl                  = 1;
        nsq_ctrls->sub_depth_block_lvl          = 2;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 10;
        nsq_ctrls->hv_weight                    = 100;
        break;
    case 29:
        nsq_ctrls->enabled                      = 1;
        nsq_ctrls->min_nsq_block_size           = 16;
        nsq_ctrls->allow_HV4                    = 0;
        nsq_ctrls->allow_HVA_HVB                = 0;
        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 5;
        nsq_ctrls->max_part0_to_part1_dev       = 0;
        nsq_ctrls->nsq_split_cost_th            = 60;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 60;
        nsq_ctrls->non_HV_split_rate_th         = 30;
        nsq_ctrls->non_HV_split_rate_modulation = 1;
        nsq_ctrls->psq_txs_lvl                  = 1;
        nsq_ctrls->sub_depth_block_lvl          = 2;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 10;
        nsq_ctrls->hv_weight                    = 100;
        break;

    case 30:
        nsq_ctrls->enabled                      = 1;
        nsq_ctrls->min_nsq_block_size           = 16;
        nsq_ctrls->allow_HV4                    = 0;
        nsq_ctrls->allow_HVA_HVB                = 0;
        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 5;
        nsq_ctrls->max_part0_to_part1_dev       = 40;
        nsq_ctrls->nsq_split_cost_th            = 60;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 60;
        nsq_ctrls->non_HV_split_rate_th         = 30;
        nsq_ctrls->non_HV_split_rate_modulation = 1;
        nsq_ctrls->psq_txs_lvl                  = 1;
        nsq_ctrls->sub_depth_block_lvl          = 2;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 10;
        nsq_ctrls->hv_weight                    = 100;
        break;

    case 31:
        nsq_ctrls->enabled                      = 1;
        nsq_ctrls->min_nsq_block_size           = 16;
        nsq_ctrls->allow_HV4                    = 0;
        nsq_ctrls->allow_HVA_HVB                = 0;
        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 5;
        nsq_ctrls->max_part0_to_part1_dev       = 60;
        nsq_ctrls->nsq_split_cost_th            = 60;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 60;
        nsq_ctrls->non_HV_split_rate_th         = 30;
        nsq_ctrls->non_HV_split_rate_modulation = 1;
        nsq_ctrls->psq_txs_lvl                  = 1;
        nsq_ctrls->sub_depth_block_lvl          = 2;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 10;
        nsq_ctrls->hv_weight                    = 100;
        break;

    case 32:
        nsq_ctrls->enabled                      = 1;
        nsq_ctrls->min_nsq_block_size           = 16;
        nsq_ctrls->allow_HV4                    = 0;
        nsq_ctrls->allow_HVA_HVB                = 0;
        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 5;
        nsq_ctrls->max_part0_to_part1_dev       = 0;
        nsq_ctrls->nsq_split_cost_th            = 60;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 60;
        nsq_ctrls->non_HV_split_rate_th         = 30;
        nsq_ctrls->non_HV_split_rate_modulation = 1;
        nsq_ctrls->psq_txs_lvl                  = 1;
        nsq_ctrls->sub_depth_block_lvl          = 2;
        nsq_ctrls->psq_pred_lvl                 = 4;
        nsq_ctrls->component_multiple_th        = 10;
        nsq_ctrls->hv_weight                    = 100;
        break;

    case 33:
        nsq_ctrls->enabled                      = 1;
        nsq_ctrls->min_nsq_block_size           = 16;
        nsq_ctrls->allow_HV4                    = 0;
        nsq_ctrls->allow_HVA_HVB                = 0;
        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 5;
        nsq_ctrls->max_part0_to_part1_dev       = 40;
        nsq_ctrls->nsq_split_cost_th            = 60;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 60;
        nsq_ctrls->non_HV_split_rate_th         = 30;
        nsq_ctrls->non_HV_split_rate_modulation = 1;
        nsq_ctrls->psq_txs_lvl                  = 1;
        nsq_ctrls->sub_depth_block_lvl          = 2;
        nsq_ctrls->psq_pred_lvl                 = 4;
        nsq_ctrls->component_multiple_th        = 10;
        nsq_ctrls->hv_weight                    = 100;
        break;

    case 34:
        nsq_ctrls->enabled                      = 1;
        nsq_ctrls->min_nsq_block_size           = 16;
        nsq_ctrls->allow_HV4                    = 0;
        nsq_ctrls->allow_HVA_HVB                = 0;
        nsq_ctrls->sq_weight                    = 95;
        nsq_ctrls->psq_cplx_lvl                 = 5;
        nsq_ctrls->max_part0_to_part1_dev       = 60;
        nsq_ctrls->nsq_split_cost_th            = 60;
        nsq_ctrls->lower_depth_split_cost_th    = 10;
        nsq_ctrls->H_vs_V_split_rate_th         = 60;
        nsq_ctrls->non_HV_split_rate_th         = 30;
        nsq_ctrls->non_HV_split_rate_modulation = 1;
        nsq_ctrls->psq_txs_lvl                  = 1;
        nsq_ctrls->sub_depth_block_lvl          = 2;
        nsq_ctrls->psq_pred_lvl                 = 4;
        nsq_ctrls->component_multiple_th        = 10;
        nsq_ctrls->hv_weight                    = 100;
        break;
    default: assert(0); break;
    }

    if (ctx && ctx->pd_pass == PD_PASS_0) {
        nsq_ctrls->sq_weight                    = (uint32_t)~0;
        nsq_ctrls->psq_cplx_lvl                 = 0;
        nsq_ctrls->max_part0_to_part1_dev       = 0;
        nsq_ctrls->nsq_split_cost_th            = 150;
        nsq_ctrls->lower_depth_split_cost_th    = 0;
        nsq_ctrls->H_vs_V_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_th         = 0;
        nsq_ctrls->non_HV_split_rate_modulation = 0;
        nsq_ctrls->psq_txs_lvl                  = 0;
        nsq_ctrls->psq_pred_lvl                 = 0;
        nsq_ctrls->component_multiple_th        = 0;
    }
    // Bypassing EncDec doesn't work if HVA_HVB_HV4 are enabled (for all bit depths; causes non-conformant bitstreams)
    if (ctx && (nsq_ctrls->allow_HV4 || nsq_ctrls->allow_HVA_HVB))
        ctx->bypass_encdec = 0;

    if (allow_HVA_HVB)
        *allow_HVA_HVB = nsq_ctrls->allow_HVA_HVB;
    if (allow_HV4)
        *allow_HV4 = nsq_ctrls->allow_HV4;
    if (min_nsq_bsize)
        *min_nsq_bsize = nsq_ctrls->min_nsq_block_size;
    if (ctx)
        memcpy(&ctx->nsq_ctrls, nsq_ctrls, sizeof(NsqCtrls));
}
void svt_aom_set_inter_intra_ctrls(ModeDecisionContext *ctx, uint8_t inter_intra_level) {
    InterIntraCompCtrls *ii_ctrls = &ctx->inter_intra_comp_ctrls;

    switch (inter_intra_level) {
    case 0:
        ii_ctrls->enabled        = 0;
        ii_ctrls->use_rd_model   = 0;
        ii_ctrls->wedge_mode_sq  = 0;
        ii_ctrls->wedge_mode_nsq = 0;
        break;
    case 1:
        ii_ctrls->enabled        = 1;
        ii_ctrls->use_rd_model   = 1;
        ii_ctrls->wedge_mode_sq  = 1;
        ii_ctrls->wedge_mode_nsq = 1;
        break;
    case 2:
        ii_ctrls->enabled        = 1;
        ii_ctrls->use_rd_model   = 0;
        ii_ctrls->wedge_mode_sq  = 0;
        ii_ctrls->wedge_mode_nsq = 2;
        break;
    default: assert(0); break;
    }
}
void svt_aom_set_depth_ctrls(PictureControlSet *pcs, ModeDecisionContext *ctx, uint8_t depth_level) {
    DepthCtrls *depth_ctrls        = &ctx->depth_ctrls;
    uint8_t     me_cplx_modulation = 0;

    switch (depth_level) {
    case 0:
        depth_ctrls->s_depth              = 0;
        depth_ctrls->e_depth              = 0;
        depth_ctrls->limit_max_min_to_pd0 = 0;
        break;
    case 1:
        depth_ctrls->s_depth              = -2;
        depth_ctrls->e_depth              = 2;
        depth_ctrls->limit_max_min_to_pd0 = 0;
        break;
    case 2:
        depth_ctrls->s_depth              = -1;
        depth_ctrls->e_depth              = 1;
        depth_ctrls->limit_max_min_to_pd0 = 0;
        break;
    case 3:
        depth_ctrls->s_depth              = -1;
        depth_ctrls->e_depth              = 1;
        depth_ctrls->limit_max_min_to_pd0 = 1;
        break;
    case 4:
        depth_ctrls->s_depth              = -1;
        depth_ctrls->e_depth              = 1;
        depth_ctrls->limit_max_min_to_pd0 = 1;
        me_cplx_modulation                = 1;
        break;
    default: assert(0); break;
    }
    SequenceControlSet *scs = pcs->scs;
    if (pcs->slice_type != I_SLICE && scs->seq_header.sb_size != BLOCK_128X128 && me_cplx_modulation) {
        if (pcs->ppcs->me_8x8_cost_variance[ctx->sb_index] <= pcs->avg_me_clpx) {
            bool is_wide_lband = (((pcs->avg_me_clpx - pcs->min_me_clpx) * 100) /
                                  MAX(pcs->max_me_clpx - pcs->min_me_clpx, 1)) > 5;

            if (is_wide_lband) {
                ctx->depth_ctrls.e_depth = 0;
            }
        } else {
            bool is_wide_rband = (((pcs->max_me_clpx - pcs->avg_me_clpx) * 100) /
                                  MAX(pcs->max_me_clpx - pcs->min_me_clpx, 1)) > 5;

            if (is_wide_rband) {
                ctx->depth_ctrls.s_depth = 0;
            }
        }
    }
}
static void set_lpd0_ctrls(ModeDecisionContext *ctx, uint8_t lpd0_lvl) {
    Lpd0Ctrls *ctrls = &ctx->lpd0_ctrls;
    // Light-PD0 only compatible with 8bit MD
    if (ctx->hbd_md) {
        ctx->lpd0_ctrls.pd0_level = REGULAR_PD0;
        return;
    }
    switch (lpd0_lvl) {
    case 0:
        ctrls->pd0_level = REGULAR_PD0; // Light-PD0 path not used
        break;
    case 1:
        ctrls->pd0_level                     = LPD0_LVL_0;
        ctrls->use_lpd0_detector[LPD0_LVL_0] = 0;
        break;
    case 2:
        ctrls->pd0_level                     = LPD0_LVL_1;
        ctrls->use_lpd0_detector[LPD0_LVL_0] = 0;
        ctrls->use_lpd0_detector[LPD0_LVL_1] = 0;
        break;
    case 3:
        ctrls->pd0_level                     = LPD0_LVL_2;
        ctrls->use_lpd0_detector[LPD0_LVL_0] = 0;
        ctrls->use_lpd0_detector[LPD0_LVL_1] = 0;
        ctrls->use_lpd0_detector[LPD0_LVL_2] = 0;
        break;
    case 4:
        ctrls->pd0_level                     = LPD0_LVL_3;
        ctrls->use_lpd0_detector[LPD0_LVL_0] = 0;
        ctrls->use_lpd0_detector[LPD0_LVL_1] = 0;
        ctrls->use_lpd0_detector[LPD0_LVL_2] = 0;
        ctrls->use_lpd0_detector[LPD0_LVL_3] = 0;
        break;
    case 5:
        ctrls->pd0_level                     = LPD0_LVL_4;
        ctrls->use_lpd0_detector[LPD0_LVL_0] = 0;
        ctrls->use_lpd0_detector[LPD0_LVL_1] = 0;
        ctrls->use_lpd0_detector[LPD0_LVL_2] = 0;
        ctrls->use_lpd0_detector[LPD0_LVL_3] = 0;

        // Set LPD0_LVL_3 controls
        ctrls->use_lpd0_detector[LPD0_LVL_4]       = 1;
        ctrls->use_ref_info[LPD0_LVL_4]            = 2;
        ctrls->me_8x8_cost_variance_th[LPD0_LVL_4] = 250000 << 1;
        ctrls->edge_dist_th[LPD0_LVL_4]            = (uint32_t)~0;
        ctrls->neigh_me_dist_shift[LPD0_LVL_4]     = (uint16_t)~0;
        break;
    case 6:
        ctrls->pd0_level                     = LPD0_LVL_4;
        ctrls->use_lpd0_detector[LPD0_LVL_0] = 0;
        ctrls->use_lpd0_detector[LPD0_LVL_1] = 0;
        ctrls->use_lpd0_detector[LPD0_LVL_2] = 0;
        ctrls->use_lpd0_detector[LPD0_LVL_3] = 0;
        // Set LPD0_LVL_3 controls
        ctrls->use_lpd0_detector[LPD0_LVL_4]       = 1;
        ctrls->use_ref_info[LPD0_LVL_4]            = 0;
        ctrls->me_8x8_cost_variance_th[LPD0_LVL_4] = 500000 << 1;
        ctrls->edge_dist_th[LPD0_LVL_4]            = (uint32_t)~0;
        ctrls->neigh_me_dist_shift[LPD0_LVL_4]     = (uint16_t)~0;
        break;
    case 7:
        ctrls->pd0_level = VERY_LIGHT_PD0;

        ctrls->use_lpd0_detector[LPD0_LVL_0] = 0;
        ctrls->use_lpd0_detector[LPD0_LVL_1] = 0;
        ctrls->use_lpd0_detector[LPD0_LVL_2] = 0;
        ctrls->use_lpd0_detector[LPD0_LVL_3] = 0;
        ctrls->use_lpd0_detector[LPD0_LVL_4] = 0;

        // Set VERY_LIGHT_PD0 controls
        ctrls->use_lpd0_detector[VERY_LIGHT_PD0]       = 1;
        ctrls->use_ref_info[VERY_LIGHT_PD0]            = 1;
        ctrls->me_8x8_cost_variance_th[VERY_LIGHT_PD0] = 250000;
        ctrls->edge_dist_th[VERY_LIGHT_PD0]            = 16384;
        ctrls->neigh_me_dist_shift[VERY_LIGHT_PD0]     = 5;
        break;
    case 8:
        ctrls->pd0_level = VERY_LIGHT_PD0;

        ctrls->use_lpd0_detector[LPD0_LVL_0] = 0;
        ctrls->use_lpd0_detector[LPD0_LVL_1] = 0;
        ctrls->use_lpd0_detector[LPD0_LVL_2] = 0;
        ctrls->use_lpd0_detector[LPD0_LVL_3] = 0;
        ctrls->use_lpd0_detector[LPD0_LVL_4] = 0;

        // Set VERY_LIGHT_PD0 controls
        ctrls->use_lpd0_detector[VERY_LIGHT_PD0]       = 1;
        ctrls->use_ref_info[VERY_LIGHT_PD0]            = 2;
        ctrls->me_8x8_cost_variance_th[VERY_LIGHT_PD0] = 250000;
        ctrls->edge_dist_th[VERY_LIGHT_PD0]            = 16384 * 5;
        ctrls->neigh_me_dist_shift[VERY_LIGHT_PD0]     = 5;
        break;
    default: assert(0); break;
    }
}
static void set_lpd1_ctrls(ModeDecisionContext *ctx, uint8_t lpd1_lvl) {
    Lpd1Ctrls *ctrls = &ctx->lpd1_ctrls;
    switch (lpd1_lvl) {
    case 0:
        ctrls->pd1_level = REGULAR_PD1; // Light-PD1 path not used
        break;
    case 1:
        ctrls->pd1_level = LPD1_LVL_0;

        // Set LPD1 level 0 controls
        ctrls->use_lpd1_detector[LPD1_LVL_0]       = 1;
        ctrls->use_ref_info[LPD1_LVL_0]            = 1;
        ctrls->cost_th_dist[LPD1_LVL_0]            = 25;
        ctrls->cost_th_rate[LPD1_LVL_0]            = 6000 + 40 * 500;
        ctrls->nz_coeff_th[LPD1_LVL_0]             = 16;
        ctrls->max_mv_length[LPD1_LVL_0]           = 300;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_0] = 250000 >> 3;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_0]   = 1024;
        ctrls->skip_pd0_me_shift[LPD1_LVL_0]       = 1;
        break;
    case 2:
        ctrls->pd1_level = LPD1_LVL_0;

        // Set LPD1 level 0 controls
        ctrls->use_lpd1_detector[LPD1_LVL_0]       = 1;
        ctrls->use_ref_info[LPD1_LVL_0]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_0]            = 35;
        ctrls->cost_th_rate[LPD1_LVL_0]            = 6000 + 100 * 500;
        ctrls->nz_coeff_th[LPD1_LVL_0]             = 16;
        ctrls->max_mv_length[LPD1_LVL_0]           = 900;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_0] = 750000 >> 3;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_0]   = 16384;
        ctrls->skip_pd0_me_shift[LPD1_LVL_0]       = 3;
        break;
    case 3:
        ctrls->pd1_level = LPD1_LVL_1;

        // Set LPD1 level 0 controls
        ctrls->use_lpd1_detector[LPD1_LVL_0]       = 1;
        ctrls->use_ref_info[LPD1_LVL_0]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_0]            = 256 << 4;
        ctrls->cost_th_rate[LPD1_LVL_0]            = 6000 + 512 * 500;
        ctrls->nz_coeff_th[LPD1_LVL_0]             = 32;
        ctrls->max_mv_length[LPD1_LVL_0]           = 2048;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_0] = 500000;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_0]   = 16384;
        ctrls->skip_pd0_me_shift[LPD1_LVL_0]       = 3;

        // Set LPD1 level 1 controls
        ctrls->use_lpd1_detector[LPD1_LVL_1]       = 1;
        ctrls->use_ref_info[LPD1_LVL_1]            = 1;
        ctrls->cost_th_dist[LPD1_LVL_1]            = 256 << 1;
        ctrls->cost_th_rate[LPD1_LVL_1]            = 6000 + 125 * 500;
        ctrls->nz_coeff_th[LPD1_LVL_1]             = 32;
        ctrls->max_mv_length[LPD1_LVL_1]           = 1600;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_1] = 500000 >> 3;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_1]   = 16384;
        ctrls->skip_pd0_me_shift[LPD1_LVL_1]       = 3;
        break;
    case 4:
        ctrls->pd1_level = LPD1_LVL_3;

        // Set LPD1 level 0 controls
        ctrls->use_lpd1_detector[LPD1_LVL_0]       = 1;
        ctrls->use_ref_info[LPD1_LVL_0]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_0]            = 256 << 10;
        ctrls->cost_th_rate[LPD1_LVL_0]            = 6000 + 8192 * 500;
        ctrls->nz_coeff_th[LPD1_LVL_0]             = 512;
        ctrls->max_mv_length[LPD1_LVL_0]           = 2048 * 16;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_0] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_0]   = 16384 * 7;
        ctrls->skip_pd0_me_shift[LPD1_LVL_0]       = 5;

        // Set LPD1 level 1 controls
        ctrls->use_lpd1_detector[LPD1_LVL_1]       = 1;
        ctrls->use_ref_info[LPD1_LVL_1]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_1]            = 256 << 8;
        ctrls->cost_th_rate[LPD1_LVL_1]            = 6000 + 4096 * 500;
        ctrls->nz_coeff_th[LPD1_LVL_1]             = 256;
        ctrls->max_mv_length[LPD1_LVL_1]           = 2048 * 8;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_1] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_1]   = 16384 * 6;
        ctrls->skip_pd0_me_shift[LPD1_LVL_1]       = 5;

        // Set LPD1 level 2 controls
        ctrls->use_lpd1_detector[LPD1_LVL_2]       = 1;
        ctrls->use_ref_info[LPD1_LVL_2]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_2]            = 256 << 8;
        ctrls->cost_th_rate[LPD1_LVL_2]            = 6000 + 4096 * 500;
        ctrls->nz_coeff_th[LPD1_LVL_2]             = 164;
        ctrls->max_mv_length[LPD1_LVL_2]           = 2048 * 8;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_2] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_2]   = 16384 * 6;
        ctrls->skip_pd0_me_shift[LPD1_LVL_2]       = 5;

        // Set LPD1 level 3 controls
        ctrls->use_lpd1_detector[LPD1_LVL_3]       = 1;
        ctrls->use_ref_info[LPD1_LVL_3]            = 1;
        ctrls->cost_th_dist[LPD1_LVL_3]            = 256 << 8;
        ctrls->cost_th_rate[LPD1_LVL_3]            = 6000 + 4096 * 500;
        ctrls->nz_coeff_th[LPD1_LVL_3]             = 128;
        ctrls->max_mv_length[LPD1_LVL_3]           = 2048 * 8;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_3] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_3]   = 16384 * 6;
        ctrls->skip_pd0_me_shift[LPD1_LVL_3]       = 5;
        break;
    case 5:
        ctrls->pd1_level = LPD1_LVL_4;

        // Set LPD1 level 0 controls
        ctrls->use_lpd1_detector[LPD1_LVL_0]       = 1;
        ctrls->use_ref_info[LPD1_LVL_0]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_0]            = 256 << 10;
        ctrls->cost_th_rate[LPD1_LVL_0]            = 6000 + 8192 * 500;
        ctrls->nz_coeff_th[LPD1_LVL_0]             = 512;
        ctrls->max_mv_length[LPD1_LVL_0]           = 2048 * 16;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_0] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_0]   = 16384 * 7;
        ctrls->skip_pd0_me_shift[LPD1_LVL_0]       = 5;

        // Set LPD1 level 1 controls
        ctrls->use_lpd1_detector[LPD1_LVL_1]       = 1;
        ctrls->use_ref_info[LPD1_LVL_1]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_1]            = 256 << 10;
        ctrls->cost_th_rate[LPD1_LVL_1]            = 6000 + 8192 * 500;
        ctrls->nz_coeff_th[LPD1_LVL_1]             = 512;
        ctrls->max_mv_length[LPD1_LVL_1]           = 2048 * 16;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_1] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_1]   = 16384 * 7;
        ctrls->skip_pd0_me_shift[LPD1_LVL_1]       = 5;

        // Set LPD1 level 2 controls
        ctrls->use_lpd1_detector[LPD1_LVL_2]       = 1;
        ctrls->use_ref_info[LPD1_LVL_2]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_2]            = 256 << 10;
        ctrls->cost_th_rate[LPD1_LVL_2]            = 6000 + 8192 * 500;
        ctrls->nz_coeff_th[LPD1_LVL_2]             = 256;
        ctrls->max_mv_length[LPD1_LVL_2]           = 2048 * 16;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_2] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_2]   = 16384 * 7;
        ctrls->skip_pd0_me_shift[LPD1_LVL_2]       = 5;

        // Set LPD1 level 3 controls
        ctrls->use_lpd1_detector[LPD1_LVL_3]       = 1;
        ctrls->use_ref_info[LPD1_LVL_3]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_3]            = 256 << 10;
        ctrls->cost_th_rate[LPD1_LVL_3]            = 6000 + 8192 * 500;
        ctrls->nz_coeff_th[LPD1_LVL_3]             = 164;
        ctrls->max_mv_length[LPD1_LVL_3]           = 2048 * 16;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_3] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_3]   = 16384 * 7;
        ctrls->skip_pd0_me_shift[LPD1_LVL_3]       = 5;

        // Set LPD1 level 4 controls
        ctrls->use_lpd1_detector[LPD1_LVL_4]       = 1;
        ctrls->use_ref_info[LPD1_LVL_4]            = 1;
        ctrls->cost_th_dist[LPD1_LVL_4]            = 256 << 4;
        ctrls->cost_th_rate[LPD1_LVL_4]            = 6000 + 1024 * 500;
        ctrls->nz_coeff_th[LPD1_LVL_4]             = 128;
        ctrls->max_mv_length[LPD1_LVL_4]           = 2048;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_4] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_4]   = 16384 * 2;
        ctrls->skip_pd0_me_shift[LPD1_LVL_4]       = 3;
        break;
    case 6:
        ctrls->pd1_level = LPD1_LVL_5;

        // Set LPD1 level 0 controls
        ctrls->use_lpd1_detector[LPD1_LVL_0]       = 1;
        ctrls->use_ref_info[LPD1_LVL_0]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_0]            = 256 << 10;
        ctrls->cost_th_rate[LPD1_LVL_0]            = 6000 + 8192 * 500;
        ctrls->nz_coeff_th[LPD1_LVL_0]             = 256;
        ctrls->max_mv_length[LPD1_LVL_0]           = 2048 * 16;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_0] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_0]   = (uint32_t)~0;
        ctrls->skip_pd0_me_shift[LPD1_LVL_0]       = (uint16_t)~0;

        // Set LPD1 level 1 controls
        ctrls->use_lpd1_detector[LPD1_LVL_1]       = 1;
        ctrls->use_ref_info[LPD1_LVL_1]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_1]            = 256 << 10;
        ctrls->cost_th_rate[LPD1_LVL_1]            = 6000 + 8192 * 500;
        ctrls->nz_coeff_th[LPD1_LVL_1]             = 128;
        ctrls->max_mv_length[LPD1_LVL_1]           = 2048 * 16;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_1] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_1]   = (uint32_t)~0;
        ctrls->skip_pd0_me_shift[LPD1_LVL_1]       = (uint16_t)~0;

        // Set LPD1 level 2 controls
        ctrls->use_lpd1_detector[LPD1_LVL_2]       = 1;
        ctrls->use_ref_info[LPD1_LVL_2]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_2]            = 256 << 10;
        ctrls->cost_th_rate[LPD1_LVL_2]            = 6000 + 8192 * 500;
        ctrls->nz_coeff_th[LPD1_LVL_2]             = 96;
        ctrls->max_mv_length[LPD1_LVL_2]           = 2048 * 16;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_2] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_2]   = (uint32_t)~0;
        ctrls->skip_pd0_me_shift[LPD1_LVL_2]       = (uint16_t)~0;

        // Set LPD1 level 3 controls
        ctrls->use_lpd1_detector[LPD1_LVL_3]       = 1;
        ctrls->use_ref_info[LPD1_LVL_3]            = 0;
        ctrls->cost_th_dist[LPD1_LVL_3]            = 256 << 10;
        ctrls->cost_th_rate[LPD1_LVL_3]            = 6000 + 8192 * 500;
        ctrls->nz_coeff_th[LPD1_LVL_3]             = 96;
        ctrls->max_mv_length[LPD1_LVL_3]           = 2048 * 16;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_3] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_3]   = (uint32_t)~0;
        ctrls->skip_pd0_me_shift[LPD1_LVL_3]       = (uint16_t)~0;

        // Set LPD1 level 4 controls
        ctrls->use_lpd1_detector[LPD1_LVL_4]       = 1;
        ctrls->use_ref_info[LPD1_LVL_4]            = 1;
        ctrls->cost_th_dist[LPD1_LVL_4]            = 256 << 8;
        ctrls->cost_th_rate[LPD1_LVL_4]            = 6000 + 4096 * 500;
        ctrls->nz_coeff_th[LPD1_LVL_4]             = 64;
        ctrls->max_mv_length[LPD1_LVL_4]           = 2048 * 8;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_4] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_4]   = 16384 * 6;
        ctrls->skip_pd0_me_shift[LPD1_LVL_4]       = 5;

        // Set LPD1 level 5 controls
        ctrls->use_lpd1_detector[LPD1_LVL_5]       = 1;
        ctrls->use_ref_info[LPD1_LVL_5]            = 1;
        ctrls->cost_th_dist[LPD1_LVL_5]            = 256 << 8;
        ctrls->cost_th_rate[LPD1_LVL_5]            = 6000 + 4096 * 500;
        ctrls->nz_coeff_th[LPD1_LVL_5]             = 64;
        ctrls->max_mv_length[LPD1_LVL_5]           = 2048 * 8;
        ctrls->me_8x8_cost_variance_th[LPD1_LVL_5] = (uint32_t)~0;
        ctrls->skip_pd0_edge_dist_th[LPD1_LVL_5]   = 16384 * 6;
        ctrls->skip_pd0_me_shift[LPD1_LVL_5]       = 5;
        break;
    default: assert(0); break;
    }
}

static void set_detect_high_freq_ctrls(ModeDecisionContext *ctx, uint8_t detect_high_freq_lvl) {
    DetectHighFreqCtrls *ctrls = &ctx->detect_high_freq_ctrls;
    switch (detect_high_freq_lvl) {
    case 0: ctrls->enabled = 0; break;
    case 1:
        ctrls->enabled             = 1;
        ctrls->high_satd_th        = 10000;
        ctrls->satd_to_sad_dev_th  = 500;
        ctrls->me_8x8_sad_var_th   = 2500;
        ctrls->depth_removal_shift = 4;
        ctrls->max_pic_lpd0_lvl    = 3;
        ctrls->max_pic_lpd1_lvl    = 2;
        ctrls->max_pd1_txt_lvl     = 8;
        break;
    case 2:
        ctrls->enabled             = 1;
        ctrls->high_satd_th        = 15000;
        ctrls->satd_to_sad_dev_th  = 600;
        ctrls->me_8x8_sad_var_th   = 7500;
        ctrls->depth_removal_shift = 4;
        ctrls->max_pic_lpd0_lvl    = 3;
        ctrls->max_pic_lpd1_lvl    = 2;
        ctrls->max_pd1_txt_lvl     = 8;
        break;
    default: assert(0); break;
    }
}
/*
 * Generate per-SB/per-PD MD settings
 */
void svt_aom_set_dist_based_ref_pruning_controls(ModeDecisionContext *ctx, uint8_t dist_based_ref_pruning_level) {
    RefPruningControls *ref_pruning_ctrls = &ctx->ref_pruning_ctrls;

    switch (dist_based_ref_pruning_level) {
    case 0: ref_pruning_ctrls->enabled = 0; break;
    case 1:
        ref_pruning_ctrls->enabled = 1;

        ref_pruning_ctrls->max_dev_to_best[PA_ME_GROUP]         = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[UNI_3x3_GROUP]       = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[BI_3x3_GROUP]        = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEW_NEAR_GROUP] = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEAR_GROUP]     = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[PRED_ME_GROUP]       = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[GLOBAL_GROUP]        = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[WARP_GROUP]          = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[OBMC_GROUP]          = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[INTER_INTRA_GROUP]   = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIST]           = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIFF]           = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[COMP_WEDGE]          = (uint32_t)~0;

        ref_pruning_ctrls->use_tpl_info_offset      = 0;
        ref_pruning_ctrls->check_closest_multiplier = 0;

        ref_pruning_ctrls->closest_refs[PA_ME_GROUP]         = 1;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP]     = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[GLOBAL_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[WARP_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[OBMC_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[INTER_INTRA_GROUP]   = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIST]           = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIFF]           = 1;
        ref_pruning_ctrls->closest_refs[COMP_WEDGE]          = 1;

        break;

    case 2:
        ref_pruning_ctrls->enabled = 1;

        ref_pruning_ctrls->max_dev_to_best[PA_ME_GROUP]         = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[UNI_3x3_GROUP]       = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[BI_3x3_GROUP]        = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEW_NEAR_GROUP] = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEAR_GROUP]     = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[PRED_ME_GROUP]       = 150;
        ref_pruning_ctrls->max_dev_to_best[GLOBAL_GROUP]        = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[WARP_GROUP]          = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[OBMC_GROUP]          = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[INTER_INTRA_GROUP]   = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIST]           = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIFF]           = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[COMP_WEDGE]          = (uint32_t)~0;

        ref_pruning_ctrls->use_tpl_info_offset      = 0;
        ref_pruning_ctrls->check_closest_multiplier = 0;

        ref_pruning_ctrls->closest_refs[PA_ME_GROUP]         = 1;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP]     = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[GLOBAL_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[WARP_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[OBMC_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[INTER_INTRA_GROUP]   = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIST]           = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIFF]           = 1;
        ref_pruning_ctrls->closest_refs[COMP_WEDGE]          = 1;

        break;
    case 3:
        ref_pruning_ctrls->enabled = 1;

        ref_pruning_ctrls->max_dev_to_best[PA_ME_GROUP]         = 30;
        ref_pruning_ctrls->max_dev_to_best[UNI_3x3_GROUP]       = 30;
        ref_pruning_ctrls->max_dev_to_best[BI_3x3_GROUP]        = 30;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEW_NEAR_GROUP] = 30;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEAR_GROUP]     = 60;
        ref_pruning_ctrls->max_dev_to_best[PRED_ME_GROUP]       = 60;
        ref_pruning_ctrls->max_dev_to_best[GLOBAL_GROUP]        = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[WARP_GROUP]          = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[OBMC_GROUP]          = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[INTER_INTRA_GROUP]   = 30;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIST]           = 30;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIFF]           = 30;
        ref_pruning_ctrls->max_dev_to_best[COMP_WEDGE]          = 30;

        ref_pruning_ctrls->use_tpl_info_offset      = 0;
        ref_pruning_ctrls->check_closest_multiplier = 0;

        ref_pruning_ctrls->closest_refs[PA_ME_GROUP]         = 1;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP]     = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[GLOBAL_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[WARP_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[OBMC_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[INTER_INTRA_GROUP]   = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIST]           = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIFF]           = 1;
        ref_pruning_ctrls->closest_refs[COMP_WEDGE]          = 1;

        break;
    case 4:
        ref_pruning_ctrls->enabled = 1;

        ref_pruning_ctrls->max_dev_to_best[PA_ME_GROUP]         = 30;
        ref_pruning_ctrls->max_dev_to_best[UNI_3x3_GROUP]       = 30;
        ref_pruning_ctrls->max_dev_to_best[BI_3x3_GROUP]        = 30;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEW_NEAR_GROUP] = 30;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEAR_GROUP]     = 30;
        ref_pruning_ctrls->max_dev_to_best[PRED_ME_GROUP]       = 30;
        ref_pruning_ctrls->max_dev_to_best[GLOBAL_GROUP]        = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[WARP_GROUP]          = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[OBMC_GROUP]          = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[INTER_INTRA_GROUP]   = 30;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIST]           = 30;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIFF]           = 30;
        ref_pruning_ctrls->max_dev_to_best[COMP_WEDGE]          = 30;

        ref_pruning_ctrls->use_tpl_info_offset      = 20;
        ref_pruning_ctrls->check_closest_multiplier = 1;

        ref_pruning_ctrls->closest_refs[PA_ME_GROUP]         = 1;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP]     = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[GLOBAL_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[WARP_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[OBMC_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[INTER_INTRA_GROUP]   = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIST]           = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIFF]           = 1;
        ref_pruning_ctrls->closest_refs[COMP_WEDGE]          = 1;

        break;
    case 5:
        ref_pruning_ctrls->enabled = 1;

        ref_pruning_ctrls->max_dev_to_best[PA_ME_GROUP]         = 30;
        ref_pruning_ctrls->max_dev_to_best[UNI_3x3_GROUP]       = 30;
        ref_pruning_ctrls->max_dev_to_best[BI_3x3_GROUP]        = 30;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEW_NEAR_GROUP] = 30;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEAR_GROUP]     = 30;
        ref_pruning_ctrls->max_dev_to_best[PRED_ME_GROUP]       = 30;
        ref_pruning_ctrls->max_dev_to_best[GLOBAL_GROUP]        = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[WARP_GROUP]          = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[OBMC_GROUP]          = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[INTER_INTRA_GROUP]   = 30;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIST]           = 0;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIFF]           = 0;
        ref_pruning_ctrls->max_dev_to_best[COMP_WEDGE]          = 0;

        ref_pruning_ctrls->use_tpl_info_offset      = 20;
        ref_pruning_ctrls->check_closest_multiplier = 1;

        ref_pruning_ctrls->closest_refs[PA_ME_GROUP]         = 1;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP]     = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[GLOBAL_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[WARP_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[OBMC_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[INTER_INTRA_GROUP]   = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIST]           = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIFF]           = 1;
        ref_pruning_ctrls->closest_refs[COMP_WEDGE]          = 1;

        break;
    case 6:
        ref_pruning_ctrls->enabled = 1;

        ref_pruning_ctrls->max_dev_to_best[PA_ME_GROUP]         = 30;
        ref_pruning_ctrls->max_dev_to_best[UNI_3x3_GROUP]       = 30;
        ref_pruning_ctrls->max_dev_to_best[BI_3x3_GROUP]        = 30;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEW_NEAR_GROUP] = 30;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEAR_GROUP]     = 30;
        ref_pruning_ctrls->max_dev_to_best[PRED_ME_GROUP]       = 30;
        ref_pruning_ctrls->max_dev_to_best[GLOBAL_GROUP]        = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[WARP_GROUP]          = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[OBMC_GROUP]          = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[INTER_INTRA_GROUP]   = 30;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIST]           = 0;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIFF]           = 0;
        ref_pruning_ctrls->max_dev_to_best[COMP_WEDGE]          = 0;

        ref_pruning_ctrls->use_tpl_info_offset      = 20;
        ref_pruning_ctrls->check_closest_multiplier = 1;

        ref_pruning_ctrls->closest_refs[PA_ME_GROUP]         = 1;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP]     = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[GLOBAL_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[WARP_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[OBMC_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[INTER_INTRA_GROUP]   = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIST]           = 0;
        ref_pruning_ctrls->closest_refs[COMP_DIFF]           = 0;
        ref_pruning_ctrls->closest_refs[COMP_WEDGE]          = 1;
        break;
    case 7:
        ref_pruning_ctrls->enabled = 1;

        ref_pruning_ctrls->max_dev_to_best[PA_ME_GROUP]         = 10;
        ref_pruning_ctrls->max_dev_to_best[UNI_3x3_GROUP]       = 10;
        ref_pruning_ctrls->max_dev_to_best[BI_3x3_GROUP]        = 10;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEW_NEAR_GROUP] = 10;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEAR_GROUP]     = 10;
        ref_pruning_ctrls->max_dev_to_best[PRED_ME_GROUP]       = 10;
        ref_pruning_ctrls->max_dev_to_best[GLOBAL_GROUP]        = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[WARP_GROUP]          = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[OBMC_GROUP]          = 10;
        ref_pruning_ctrls->max_dev_to_best[INTER_INTRA_GROUP]   = 10;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIST]           = 0;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIFF]           = 0;
        ref_pruning_ctrls->max_dev_to_best[COMP_WEDGE]          = 0;

        ref_pruning_ctrls->use_tpl_info_offset      = 20;
        ref_pruning_ctrls->check_closest_multiplier = 1;

        ref_pruning_ctrls->closest_refs[PA_ME_GROUP]         = 1;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP]     = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[GLOBAL_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[WARP_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[OBMC_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[INTER_INTRA_GROUP]   = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIST]           = 0;
        ref_pruning_ctrls->closest_refs[COMP_DIFF]           = 0;
        ref_pruning_ctrls->closest_refs[COMP_WEDGE]          = 1;

        break;
    case 8:
        ref_pruning_ctrls->enabled = 1;

        ref_pruning_ctrls->max_dev_to_best[PA_ME_GROUP]         = 0;
        ref_pruning_ctrls->max_dev_to_best[UNI_3x3_GROUP]       = 0;
        ref_pruning_ctrls->max_dev_to_best[BI_3x3_GROUP]        = 0;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEW_NEAR_GROUP] = 0;
        ref_pruning_ctrls->max_dev_to_best[NRST_NEAR_GROUP]     = 0;
        ref_pruning_ctrls->max_dev_to_best[PRED_ME_GROUP]       = 0;
        ref_pruning_ctrls->max_dev_to_best[GLOBAL_GROUP]        = (uint32_t)~0;
        ref_pruning_ctrls->max_dev_to_best[WARP_GROUP]          = 0;
        ref_pruning_ctrls->max_dev_to_best[OBMC_GROUP]          = 0;
        ref_pruning_ctrls->max_dev_to_best[INTER_INTRA_GROUP]   = 0;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIST]           = 0;
        ref_pruning_ctrls->max_dev_to_best[COMP_DIFF]           = 0;
        ref_pruning_ctrls->max_dev_to_best[COMP_WEDGE]          = 0;

        ref_pruning_ctrls->use_tpl_info_offset      = 20;
        ref_pruning_ctrls->check_closest_multiplier = 1;

        ref_pruning_ctrls->closest_refs[PA_ME_GROUP]         = 1;
        ref_pruning_ctrls->closest_refs[UNI_3x3_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[BI_3x3_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEW_NEAR_GROUP] = 1;
        ref_pruning_ctrls->closest_refs[NRST_NEAR_GROUP]     = 1;
        ref_pruning_ctrls->closest_refs[PRED_ME_GROUP]       = 1;
        ref_pruning_ctrls->closest_refs[GLOBAL_GROUP]        = 1;
        ref_pruning_ctrls->closest_refs[WARP_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[OBMC_GROUP]          = 1;
        ref_pruning_ctrls->closest_refs[INTER_INTRA_GROUP]   = 1;
        ref_pruning_ctrls->closest_refs[COMP_DIST]           = 0;
        ref_pruning_ctrls->closest_refs[COMP_DIFF]           = 0;
        ref_pruning_ctrls->closest_refs[COMP_WEDGE]          = 1;
        break;
    default: assert(0); break;
    }
}
static void set_txs_controls(ModeDecisionContext *ctx, uint8_t txs_level) {
    TxsControls *txs_ctrls = &ctx->txs_ctrls;

    switch (txs_level) {
    case 0: txs_ctrls->enabled = 0; break;
    case 1:
        txs_ctrls->enabled                   = 1;
        txs_ctrls->prev_depth_coeff_exit_th  = 1;
        txs_ctrls->intra_class_max_depth_sq  = 2;
        txs_ctrls->intra_class_max_depth_nsq = 2;
        txs_ctrls->inter_class_max_depth_sq  = 2;
        txs_ctrls->inter_class_max_depth_nsq = 2;
        txs_ctrls->depth1_txt_group_offset   = 0;
        txs_ctrls->depth2_txt_group_offset   = 0;
        txs_ctrls->min_sq_size               = 0;
        break;
    case 2:
        txs_ctrls->enabled                   = 1;
        txs_ctrls->prev_depth_coeff_exit_th  = 1;
        txs_ctrls->intra_class_max_depth_sq  = 2;
        txs_ctrls->intra_class_max_depth_nsq = 2;
        txs_ctrls->inter_class_max_depth_sq  = 1;
        txs_ctrls->inter_class_max_depth_nsq = 1;
        txs_ctrls->depth1_txt_group_offset   = 0;
        txs_ctrls->depth2_txt_group_offset   = 0;
        txs_ctrls->min_sq_size               = 0;
        break;
    case 3:
        txs_ctrls->enabled                   = 1;
        txs_ctrls->prev_depth_coeff_exit_th  = 1;
        txs_ctrls->intra_class_max_depth_sq  = 2;
        txs_ctrls->intra_class_max_depth_nsq = 1;
        txs_ctrls->inter_class_max_depth_sq  = 1;
        txs_ctrls->inter_class_max_depth_nsq = 0;
        txs_ctrls->depth1_txt_group_offset   = 0;
        txs_ctrls->depth2_txt_group_offset   = 0;
        txs_ctrls->min_sq_size               = 0;
        break;
    case 4:
        txs_ctrls->enabled                   = 1;
        txs_ctrls->prev_depth_coeff_exit_th  = 1;
        txs_ctrls->intra_class_max_depth_sq  = 1;
        txs_ctrls->intra_class_max_depth_nsq = 0;
        txs_ctrls->inter_class_max_depth_sq  = 1;
        txs_ctrls->inter_class_max_depth_nsq = 0;
        txs_ctrls->depth1_txt_group_offset   = 4;
        txs_ctrls->depth2_txt_group_offset   = 4;
        txs_ctrls->min_sq_size               = 0;
        break;
    case 5:
        txs_ctrls->enabled                   = 1;
        txs_ctrls->prev_depth_coeff_exit_th  = 1;
        txs_ctrls->intra_class_max_depth_sq  = 1;
        txs_ctrls->intra_class_max_depth_nsq = 0;
        txs_ctrls->inter_class_max_depth_sq  = 0;
        txs_ctrls->inter_class_max_depth_nsq = 0;
        txs_ctrls->depth1_txt_group_offset   = 4;
        txs_ctrls->depth2_txt_group_offset   = 4;
        txs_ctrls->min_sq_size               = 0;
        break;
    default: assert(0); break;
    }
}
static void set_spatial_sse_full_loop_level(ModeDecisionContext *ctx, uint8_t spatial_sse_full_loop_level) {
    SpatialSSECtrls *spatial_sse_ctrls = &ctx->spatial_sse_ctrls;

    switch (spatial_sse_full_loop_level) {
    case 0: spatial_sse_ctrls->spatial_sse_full_loop_level = FALSE; break;
    case 1: spatial_sse_ctrls->spatial_sse_full_loop_level = TRUE; break;
    default: assert(0); break;
    }
}
// Compute a qp-aware threshold based on the variance of the SB, used to apply selectively INTRA at PD0
static uint64_t compute_intra_pd0_th(SequenceControlSet *scs, PictureControlSet *pcs, ModeDecisionContext *ctx) {
    uint32_t fast_lambda      = ctx->hbd_md ? ctx->fast_lambda_md[EB_10_BIT_MD] : ctx->fast_lambda_md[EB_8_BIT_MD];
    uint32_t sb_size          = scs->super_block_size * scs->super_block_size;
    uint64_t cost_th_rate     = 1 << 13;
    uint64_t use_intra_pd0_th = 0;

    if (scs->calculate_variance) {
        const uint16_t var_th0 = 400; // low-texture 64x64
        const uint16_t var_th1 = 800; // moderate-texture 64x64
        if (pcs->ppcs->variance[ctx->sb_index][ME_TIER_ZERO_PU_64x64] <= var_th0)
            use_intra_pd0_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 8);
        else if (pcs->ppcs->variance[ctx->sb_index][ME_TIER_ZERO_PU_64x64] <= var_th1)
            use_intra_pd0_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 7);
        else
            use_intra_pd0_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 6);
    } else {
        use_intra_pd0_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 6);
    }
    return use_intra_pd0_th;
}
// Compute a qp-aware threshold based on the variance of the SB, used to apply selectively subres
static uint64_t compute_subres_th(SequenceControlSet *scs, PictureControlSet *pcs, ModeDecisionContext *ctx) {
    uint32_t fast_lambda   = ctx->hbd_md ? ctx->fast_lambda_md[EB_10_BIT_MD] : ctx->fast_lambda_md[EB_8_BIT_MD];
    uint32_t sb_size       = scs->super_block_size * scs->super_block_size;
    uint64_t cost_th_rate  = 1 << 13;
    uint64_t use_subres_th = 0;

    if (scs->calculate_variance) {
        if (pcs->ppcs->variance[ctx->sb_index][ME_TIER_ZERO_PU_64x64] <= 400)
            use_subres_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 8);
        else if (pcs->ppcs->variance[ctx->sb_index][ME_TIER_ZERO_PU_64x64] <= 800)
            use_subres_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 7);
        else
            use_subres_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 6);
    } else {
        use_subres_th = RDCOST(fast_lambda, cost_th_rate, sb_size * 6);
    }
    return use_subres_th;
}
static void set_lpd1_tx_ctrls(PictureControlSet *pcs, ModeDecisionContext *ctx, uint8_t lpd1_tx_level) {
    PictureParentControlSet *ppcs  = pcs->ppcs;
    Lpd1TxCtrls             *ctrls = &ctx->lpd1_tx_ctrls;

    switch (lpd1_tx_level) {
    case 0:
        ctrls->zero_y_coeff_exit            = 0;
        ctrls->skip_nrst_nrst_luma_tx       = 0;
        ctrls->skip_tx_th                   = 0;
        ctrls->use_uv_shortcuts_on_y_coeffs = 0;

        ctrls->use_mds3_shortcuts_th = 0;
        ctrls->use_neighbour_info    = 0;
        break;
    case 1:
        ctrls->zero_y_coeff_exit            = 1;
        ctrls->chroma_detector_level        = 1;
        ctrls->skip_nrst_nrst_luma_tx       = 0;
        ctrls->skip_tx_th                   = 0;
        ctrls->use_uv_shortcuts_on_y_coeffs = 1;

        ctrls->use_mds3_shortcuts_th = 25;
        ctrls->use_neighbour_info    = 0;
        break;
    case 2:
        ctrls->zero_y_coeff_exit            = 1;
        ctrls->chroma_detector_level        = 1;
        ctrls->skip_nrst_nrst_luma_tx       = 1;
        ctrls->skip_tx_th                   = 25;
        ctrls->use_uv_shortcuts_on_y_coeffs = 1;

        ctrls->use_mds3_shortcuts_th = 25;
        ctrls->use_neighbour_info    = 0;
        break;
    case 3:
        ctrls->zero_y_coeff_exit            = 1;
        ctrls->chroma_detector_level        = ppcs->is_ref ? 2 : 3;
        ctrls->skip_nrst_nrst_luma_tx       = 1;
        ctrls->skip_tx_th                   = 25;
        ctrls->use_uv_shortcuts_on_y_coeffs = 1;

        ctrls->use_mds3_shortcuts_th = 25;
        ctrls->use_neighbour_info    = 0;
        break;
    case 4:
        ctrls->zero_y_coeff_exit            = 1;
        ctrls->chroma_detector_level        = ppcs->is_ref ? 2 : 3;
        ctrls->skip_nrst_nrst_luma_tx       = 1;
        ctrls->skip_tx_th                   = 25;
        ctrls->use_uv_shortcuts_on_y_coeffs = 1;

        ctrls->use_mds3_shortcuts_th = 50;
        ctrls->use_neighbour_info    = 1;
        break;
    case 5:
        ctrls->zero_y_coeff_exit            = 1;
        ctrls->chroma_detector_level        = ppcs->is_ref ? 2 : 3;
        ctrls->skip_nrst_nrst_luma_tx       = 1;
        ctrls->skip_tx_th                   = 50;
        ctrls->use_uv_shortcuts_on_y_coeffs = 1;

        ctrls->use_mds3_shortcuts_th = 50;
        ctrls->use_neighbour_info    = 2;
        break;
    case 6:
        ctrls->zero_y_coeff_exit            = 1;
        ctrls->chroma_detector_level        = ppcs->is_ref ? 2 : 3;
        ctrls->skip_nrst_nrst_luma_tx       = 1;
        ctrls->skip_tx_th                   = 70;
        ctrls->use_uv_shortcuts_on_y_coeffs = 1;

        ctrls->use_mds3_shortcuts_th = 100;
        ctrls->use_neighbour_info    = 2;
        break;
    default: assert(0); break;
    }
}
static void set_cfl_ctrls(ModeDecisionContext *ctx, uint8_t cfl_level) {
    CflCtrls *ctrls = &ctx->cfl_ctrls;

    switch (cfl_level) {
    case 0: ctrls->enabled = 0; break;
    case 1:
        ctrls->enabled = 1;
        ctrls->itr_th  = 2;
        break;
    case 2:
        ctrls->enabled = 1;
        ctrls->itr_th  = 1;
        break;
    default: assert(0); break;
    }
}
static void set_rate_est_ctrls(ModeDecisionContext *ctx, uint8_t rate_est_level) {
    MdRateEstCtrls *ctrls = &ctx->rate_est_ctrls;

    switch (rate_est_level) {
    case 0:
        ctrls->update_skip_ctx_dc_sign_ctx = 0;
        ctrls->update_skip_coeff_ctx       = 0;
        ctrls->coeff_rate_est_lvl          = 0;
        ctrls->lpd0_qp_offset              = 8;
        ctrls->pd0_fast_coeff_est_level    = 2;
        break;
    case 1:
        ctrls->update_skip_ctx_dc_sign_ctx = 1;
        ctrls->update_skip_coeff_ctx       = 1;
        ctrls->coeff_rate_est_lvl          = 1;
        ctrls->lpd0_qp_offset              = 0;
        ctrls->pd0_fast_coeff_est_level    = 1;
        break;
    case 2:
        ctrls->update_skip_ctx_dc_sign_ctx = 1;
        ctrls->update_skip_coeff_ctx       = 0;
        ctrls->coeff_rate_est_lvl          = 1;
        ctrls->lpd0_qp_offset              = 0;
        ctrls->pd0_fast_coeff_est_level    = 2;
        break;
    case 3:
        ctrls->update_skip_ctx_dc_sign_ctx = 1;
        ctrls->update_skip_coeff_ctx       = 0;
        ctrls->coeff_rate_est_lvl          = 2;
        ctrls->lpd0_qp_offset              = 0;
        ctrls->pd0_fast_coeff_est_level    = 2;
        break;
    case 4:
        ctrls->update_skip_ctx_dc_sign_ctx = 0;
        ctrls->update_skip_coeff_ctx       = 0;
        ctrls->coeff_rate_est_lvl          = 2;
        ctrls->lpd0_qp_offset              = 0;
        ctrls->pd0_fast_coeff_est_level    = 2;
        break;
    default: assert(0); break;
    }
}
/*
Loop over TPL blocks in the SB to update intra information.  Return 1 if the stats for the SB are valid; else return 0.

sb_ang_intra_count: Number of TPL blocks in the SB where the best_mode was an angular intra mode
sb_max_intra: The maximum intra mode selected by any TPL block in the SB (DC_PRED is lowest, PAETH_PRED is highest)
sb_intra_count: Number of TPL blocks in the SB where the best_mode was an intra mode
*/
static Bool get_sb_tpl_intra_stats(PictureControlSet *pcs, ModeDecisionContext *ctx, int *sb_ang_intra_count,
                                   PredictionMode *sb_max_intra, int *sb_intra_count) {
    PictureParentControlSet *ppcs = pcs->ppcs;

    // Check that TPL data is available and that INTRA was tested in TPL.
    // Note that not all INTRA modes may be tested in TPL.
    if (ppcs->tpl_ctrls.enable && ppcs->tpl_src_data_ready &&
        (ppcs->is_ref || !ppcs->tpl_ctrls.disable_intra_pred_nref)) {
        const int      aligned16_width = (ppcs->aligned_width + 15) >> 4;
        const uint32_t mb_origin_x     = ctx->sb_origin_x;
        const uint32_t mb_origin_y     = ctx->sb_origin_y;
        const int      tpl_blk_size    = ppcs->tpl_ctrls.dispenser_search_level == 0 ? 16
                    : ppcs->tpl_ctrls.dispenser_search_level == 1                    ? 32
                                                                                     : 64;

        // Get actual SB width (for cases of incomplete SBs)
        SbGeom *sb_geom = &ppcs->sb_geom[ctx->sb_index];
        int     sb_cols = MAX(1, sb_geom->width / tpl_blk_size);
        int     sb_rows = MAX(1, sb_geom->height / tpl_blk_size);

        int            ang_intra_count = 0;
        PredictionMode max_intra       = DC_PRED;
        int            intra_count     = 0;

        // Loop over all blocks in the SB
        for (int i = 0; i < sb_rows; i++) {
            TplSrcStats *tpl_src_stats_buffer =
                &ppcs->pa_me_data
                     ->tpl_src_stats_buffer[((mb_origin_y >> 4) + i) * aligned16_width + (mb_origin_x >> 4)];
            for (int j = 0; j < sb_cols; j++) {
                if (is_intra_mode(tpl_src_stats_buffer->best_mode)) {
                    max_intra = MAX(max_intra, tpl_src_stats_buffer->best_mode);
                    intra_count++;
                }

                if (av1_is_directional_mode(tpl_src_stats_buffer->best_mode)) {
                    ang_intra_count++;
                }
                tpl_src_stats_buffer++;
            }
        }

        *sb_ang_intra_count = ang_intra_count;
        *sb_max_intra       = max_intra;
        *sb_intra_count     = intra_count;
        return 1;
    }
    return 0;
}
static void set_intra_ctrls(PictureControlSet *pcs, ModeDecisionContext *ctx, uint8_t intra_level) {
    IntraCtrls              *ctrls = &ctx->intra_ctrls;
    PictureParentControlSet *ppcs  = pcs->ppcs;

    // If intra is disallowed at the pic level, must disallow at SB level
    if (pcs->skip_intra)
        intra_level = 0;

    assert(IMPLIES(pcs->slice_type == I_SLICE, intra_level > 0));

    switch (intra_level) {
    case 0:
        ctrls->enable_intra       = 0;
        ctrls->intra_mode_end     = DC_PRED;
        ctrls->angular_pred_level = 0;
        break;
    case 1:
        ctrls->enable_intra       = 1;
        ctrls->intra_mode_end     = PAETH_PRED;
        ctrls->angular_pred_level = 1;
        break;
    case 2:
        ctrls->enable_intra       = 1;
        ctrls->intra_mode_end     = PAETH_PRED;
        ctrls->angular_pred_level = 2;

        // Only use TPL info if all INTRA modes are tested
        if (ppcs->tpl_ctrls.enable && ppcs->tpl_ctrls.intra_mode_end == PAETH_PRED) {
            int            sb_ang_intra_count;
            PredictionMode sb_max_intra;
            int            sb_intra_count;
            if (get_sb_tpl_intra_stats(pcs, ctx, &sb_ang_intra_count, &sb_max_intra, &sb_intra_count)) {
                // if SB has angluar modes, use full search
                if (sb_ang_intra_count) {
                    ctrls->angular_pred_level = 1;
                } else {
                    ctrls->angular_pred_level = 3;
                }
            }
        }
        break;
    case 3:
        ctrls->enable_intra       = 1;
        ctrls->intra_mode_end     = PAETH_PRED;
        ctrls->angular_pred_level = 2;

        // Only use TPL info if all INTRA modes are tested
        if (ppcs->tpl_ctrls.enable && ppcs->tpl_ctrls.intra_mode_end == PAETH_PRED) {
            int            sb_ang_intra_count;
            PredictionMode sb_max_intra;
            int            sb_intra_count;
            if (get_sb_tpl_intra_stats(pcs, ctx, &sb_ang_intra_count, &sb_max_intra, &sb_intra_count)) {
                // if SB has angluar modes, use full search
                if (sb_ang_intra_count) {
                    ctrls->angular_pred_level = 1;
                } else {
                    ctrls->angular_pred_level = 3;
                }
                ctrls->intra_mode_end = sb_max_intra;
            }
        }
        break;
    case 4:
        ctrls->enable_intra       = 1;
        ctrls->intra_mode_end     = SMOOTH_H_PRED;
        ctrls->angular_pred_level = 3;

        // Only use TPL info if all INTRA modes are tested
        if (ppcs->tpl_ctrls.enable && ppcs->tpl_ctrls.intra_mode_end == PAETH_PRED) {
            int            sb_ang_intra_count;
            PredictionMode sb_max_intra;
            int            sb_intra_count;
            if (get_sb_tpl_intra_stats(pcs, ctx, &sb_ang_intra_count, &sb_max_intra, &sb_intra_count)) {
                int tpl_blk_size = ppcs->tpl_ctrls.dispenser_search_level == 0 ? 16
                    : ppcs->tpl_ctrls.dispenser_search_level == 1              ? 32
                                                                               : 64;
                // Get actual SB width (for cases of incomplete SBs)
                SbGeom *sb_geom = &ppcs->sb_geom[ctx->sb_index];
                int     sb_cols = sb_geom->width / tpl_blk_size;
                int     sb_rows = sb_geom->height / tpl_blk_size;

                // if more than a quarter of SB is angular, use safe intra_level
                if (sb_ang_intra_count > ((sb_rows * sb_cols) >> 2)) {
                    ctrls->angular_pred_level = 1;
                } else if (sb_ang_intra_count > 2) {
                    ctrls->angular_pred_level = 2;
                } else {
                    ctrls->angular_pred_level = 4;
                }

                ctrls->intra_mode_end = sb_max_intra;
            }
        }
        break;
    case 5:
        ctrls->enable_intra       = 1;
        ctrls->intra_mode_end     = SMOOTH_PRED;
        ctrls->angular_pred_level = 4;

        // There is no check that all TPL modes are checked, so should only use info about
        // general intra modes, not the specific intra mode selected or whether it's angular
        if (ppcs->tpl_ctrls.enable) {
            int            sb_ang_intra_count;
            PredictionMode sb_max_intra;
            int            sb_intra_count;
            if (get_sb_tpl_intra_stats(pcs, ctx, &sb_ang_intra_count, &sb_max_intra, &sb_intra_count)) {
                if (sb_intra_count > 0) {
                    ctrls->angular_pred_level = 2;
                    ctrls->intra_mode_end     = PAETH_PRED;
                } else if (pcs->ref_intra_percentage < 30) {
                    ctrls->angular_pred_level = 0;
                    ctrls->intra_mode_end     = SMOOTH_PRED;
                }
            }
        }
        break;
    case 6:
        ctrls->enable_intra       = 1;
        ctrls->intra_mode_end     = SMOOTH_PRED;
        ctrls->angular_pred_level = 4;
        break;
    case 7:
        ctrls->enable_intra       = 1;
        ctrls->intra_mode_end     = DC_PRED;
        ctrls->angular_pred_level = 0;
        break;
    default: assert(0); break;
    }

    /* For PD1, the ability to skip intra must be set at the pic level to ensure all SBs
    perform inverse TX and generate the recon. */
    if (ctx->pd_pass == PD_PASS_1) {
        ctx->skip_intra = pcs->skip_intra;

        // Check user-defined settings
        if (pcs->ppcs->scs->enable_paeth == 0)
            ctrls->intra_mode_end = MIN(ctrls->intra_mode_end, SMOOTH_H_PRED);

        if (pcs->ppcs->scs->enable_smooth == 0)
            ctrls->intra_mode_end = MIN(ctrls->intra_mode_end, D67_PRED);

        if (pcs->ppcs->scs->intra_angle_delta == 0)
            ctrls->angular_pred_level = 0;
    } else {
        ctx->skip_intra = !(ctrls->enable_intra) || pcs->skip_intra;
        assert(IMPLIES(ctx->lpd0_ctrls.pd0_level == VERY_LIGHT_PD0, ctx->skip_intra == 1));
    }
}

static void set_tx_shortcut_ctrls(PictureControlSet *pcs, ModeDecisionContext *ctx, uint8_t tx_shortcut_level) {
    PictureParentControlSet *ppcs  = pcs->ppcs;
    TxShortcutCtrls         *ctrls = &ctx->tx_shortcut_ctrls;

    switch (tx_shortcut_level) {
    case 0:
        ctrls->bypass_tx_when_zcoeff = 0;
        ctrls->apply_pf_on_coeffs    = 0;
        ctrls->use_mds3_shortcuts_th = 0;
        ctrls->use_neighbour_info    = 0;
        ctrls->chroma_detector_level = 0;
        break;
    case 1:
        ctrls->bypass_tx_when_zcoeff = 1;
        ctrls->apply_pf_on_coeffs    = 0;
        ctrls->use_mds3_shortcuts_th = 0;
        ctrls->chroma_detector_level = 1;
        ctrls->use_neighbour_info    = 0;
        break;
    case 2:
        ctrls->bypass_tx_when_zcoeff = 1;
        ctrls->apply_pf_on_coeffs    = 1;
        ctrls->use_mds3_shortcuts_th = 0;
        ctrls->chroma_detector_level = 1;
        ctrls->use_neighbour_info    = 0;
        break;
    case 3:
        ctrls->bypass_tx_when_zcoeff = 1;
        ctrls->apply_pf_on_coeffs    = 1;
        ctrls->use_mds3_shortcuts_th = 25;
        ctrls->chroma_detector_level = 1;
        ctrls->use_neighbour_info    = 0;
        break;
    case 4:
        ctrls->bypass_tx_when_zcoeff = 1;
        ctrls->apply_pf_on_coeffs    = 1;
        ctrls->use_mds3_shortcuts_th = 25;
        ctrls->chroma_detector_level = ppcs->is_ref ? 1 : 0;
        ctrls->use_neighbour_info    = 1;
        break;
    case 5:
        ctrls->bypass_tx_when_zcoeff = 1;
        ctrls->apply_pf_on_coeffs    = 1;
        ctrls->use_mds3_shortcuts_th = 50;
        ctrls->chroma_detector_level = 0;
        ctrls->use_neighbour_info    = 1;
        break;
    default: assert(0); break;
    }

    // Chroma detector should be used in M12 and below (at least in REF frames) to prevent blurring artifacts in some clips
    if (tx_shortcut_level && ppcs->is_ref && pcs->enc_mode <= ENC_M12)
        assert(ctrls->chroma_detector_level &&
               "Chroma detector should be used for ref frames in low presets to prevent blurring "
               "artifacts.");
}

static void set_mds0_controls(ModeDecisionContext *ctx, uint8_t mds0_level) {
    Mds0Ctrls *ctrls = &ctx->mds0_ctrls;

    switch (mds0_level) {
    case 0:
        ctrls->mds0_dist_type               = MDS0_SAD;
        ctrls->enable_cost_based_early_exit = 1;
        ctrls->mds0_distortion_th           = 0;
        break;
    case 1:
        ctrls->mds0_dist_type               = MDS0_SSD;
        ctrls->enable_cost_based_early_exit = 0;
        ctrls->mds0_distortion_th           = 0;
        break;
    case 2:
        ctrls->mds0_dist_type               = MDS0_VAR;
        ctrls->enable_cost_based_early_exit = 0;
        ctrls->mds0_distortion_th           = 0;
        break;
    case 3:
        ctrls->mds0_dist_type               = MDS0_VAR;
        ctrls->enable_cost_based_early_exit = 1;
        ctrls->mds0_distortion_th           = 50;
        break;
    case 4:
        ctrls->mds0_dist_type               = MDS0_VAR;
        ctrls->enable_cost_based_early_exit = 1;
        ctrls->mds0_distortion_th           = 0;
        break;
    default: assert(0); break;
    }
}
static void set_skip_sub_depth_ctrls(SkipSubDepthCtrls *skip_sub_depth_ctrls, uint8_t skip_sub_depth_lvl) {
    switch (skip_sub_depth_lvl) {
    case 0: skip_sub_depth_ctrls->enabled = 0; break;

    case 1:
        skip_sub_depth_ctrls->enabled = 1;
        // Cond0 ctrls
        skip_sub_depth_ctrls->max_size_cond0 = 8;
        skip_sub_depth_ctrls->nsq_to_sq_th   = 5;
        // Cond1 ctrls
        skip_sub_depth_ctrls->max_size_cond1    = 0;
        skip_sub_depth_ctrls->quad_deviation_th = 250;
        skip_sub_depth_ctrls->coeff_perc        = 25;
        break;

    case 2:

        skip_sub_depth_ctrls->enabled = 1;
        // Cond0 ctrls
        skip_sub_depth_ctrls->max_size_cond0 = 8;
        skip_sub_depth_ctrls->nsq_to_sq_th   = 1;
        // Cond1 ctrls
        skip_sub_depth_ctrls->max_size_cond1    = 16;
        skip_sub_depth_ctrls->quad_deviation_th = 250;
        skip_sub_depth_ctrls->coeff_perc        = 15;

        break;

    case 3:

        skip_sub_depth_ctrls->enabled = 1;
        // Cond0 ctrls
        skip_sub_depth_ctrls->max_size_cond0 = 8;
        skip_sub_depth_ctrls->nsq_to_sq_th   = 1;
        // Cond1 ctrls
        skip_sub_depth_ctrls->max_size_cond1    = 16;
        skip_sub_depth_ctrls->quad_deviation_th = 250;
        skip_sub_depth_ctrls->coeff_perc        = 25;

        break;

    case 4:

        skip_sub_depth_ctrls->enabled = 1;
        // Cond0 ctrls
        skip_sub_depth_ctrls->max_size_cond0 = 8;
        skip_sub_depth_ctrls->nsq_to_sq_th   = 0;
        // Cond1 ctrls
        skip_sub_depth_ctrls->max_size_cond1    = 128;
        skip_sub_depth_ctrls->quad_deviation_th = 250;
        skip_sub_depth_ctrls->coeff_perc        = 25;

        break;

    case 5:

        skip_sub_depth_ctrls->enabled = 1;
        // Cond0 ctrls
        skip_sub_depth_ctrls->max_size_cond0 = 8;
        skip_sub_depth_ctrls->nsq_to_sq_th   = 0;
        // Cond1 ctrls
        skip_sub_depth_ctrls->max_size_cond1    = 128;
        skip_sub_depth_ctrls->quad_deviation_th = 250;
        skip_sub_depth_ctrls->coeff_perc        = 50;

        break;

    default: assert(0); break;
    }
}
/*
 * Generate per-SB MD settings (do not change per-PD)
 */
void svt_aom_sig_deriv_enc_dec_common(SequenceControlSet *scs, PictureControlSet *pcs, ModeDecisionContext *ctx) {
    EncMode enc_mode = pcs->enc_mode;

    SuperBlock *sb_ptr = pcs->sb_ptr_array[ctx->sb_index];
    set_detect_high_freq_ctrls(ctx, pcs->vq_ctrls.detect_high_freq_lvl);
    ctx->high_freq_present = 0;
    const bool rtc_tune    = (scs->static_config.pred_structure == SVT_AV1_PRED_LOW_DELAY_B) ? true : false;
    if (ctx->detect_high_freq_ctrls.enabled) {
        svt_aom_check_high_freq(pcs, sb_ptr, ctx);
    }

    // Level 0: pred depth only
    // Level 1: [-2, +2] depth refinement
    // Level 2: [-1, +1] depth refinement
    uint8_t depth_level = 0;
    if (pcs->ppcs->sc_class1) {
        if (enc_mode <= ENC_M4)
            depth_level = pcs->slice_type == I_SLICE ? 1 : 2;
        else if (enc_mode <= ENC_M8)
            depth_level = 2;
        else
            depth_level = pcs->slice_type == I_SLICE ? 2 : 0;
    } else if (rtc_tune) {
        if (enc_mode <= ENC_M8)
            depth_level = 2;
        else if (enc_mode <= ENC_M9) {
            depth_level = 2;
            if (pcs->slice_type != I_SLICE) {
                int me_8x8 = pcs->ppcs->me_8x8_cost_variance[ctx->sb_index];
                int th     = enc_mode <= ENC_M9 ? 60 * ctx->qp_index : 500 * ctx->qp_index;
                if (me_8x8 < th)
                    depth_level = 0;
            }
        } else if (enc_mode <= ENC_M10)
            depth_level = pcs->slice_type == I_SLICE ? 2 : 0;
        else
            depth_level = 0;
    } else if (enc_mode <= ENC_M5)
        depth_level = pcs->slice_type == I_SLICE ? 1 : 2;
    else if (enc_mode <= ENC_M7) {
        if (pcs->coeff_lvl == LOW_LVL) {
            depth_level = pcs->slice_type == I_SLICE ? 1 : 2;
        } else if (pcs->coeff_lvl == HIGH_LVL) {
            depth_level = 4;
        } else {
            depth_level = 3;
        }
    } else if (enc_mode <= ENC_M9) {
        if (pcs->coeff_lvl == LOW_LVL) {
            depth_level = 3;
        } else if (pcs->coeff_lvl == HIGH_LVL) {
            depth_level = 0;
        } else {
            depth_level = 4;
        }
    } else
        depth_level = 0;
    svt_aom_set_depth_ctrls(pcs, ctx, depth_level);
    // pic_pred_depth_only shouldn't be changed after this point
    ctx->pred_depth_only = ctx->pic_pred_depth_only = (depth_level == 0);

    set_lpd0_ctrls(ctx,
                   ctx->high_freq_present ? MIN(pcs->pic_lpd0_lvl, ctx->detect_high_freq_ctrls.max_pic_lpd0_lvl)
                                          : pcs->pic_lpd0_lvl);

    B64Geom *b64_geom                             = &pcs->ppcs->b64_geom[ctx->sb_index];
    ctx->depth_removal_ctrls.disallow_below_64x64 = 0;
    ctx->depth_removal_ctrls.disallow_below_32x32 = 0;
    /*
if disallow_below_16x16 is turned ON then enable_me_8x8 should be turned OFF for the same preset in order to save memory and cycles as that feature optimizes the me_candidate_array,
me_mv_array and the total_me_candidate_index arrays when 8x8 blocks are not used

if any check other than an I-SLICE check is used on disallow_below_16x16 then the enable_me_8x8 should be turned ON for the entire preset because without the 8x8 me data the non I-SLICE pictures
that use 8x8 blocks will lose significant BD-Rate as the parent 16x16 me data will be used for the 8x8 blocks
*/
    ctx->depth_removal_ctrls.disallow_below_16x16 = pcs->pic_disallow_below_16x16;

    if (b64_geom->width % 32 != 0 || b64_geom->height % 32 != 0)
        ctx->depth_removal_ctrls.disallow_below_64x64 = FALSE;
    if (b64_geom->width % 16 != 0 || b64_geom->height % 16 != 0)
        ctx->depth_removal_ctrls.disallow_below_32x32 = FALSE;
    if (b64_geom->width % 8 != 0 || b64_geom->height % 8 != 0)
        ctx->depth_removal_ctrls.disallow_below_16x16 = FALSE;
    ctx->disallow_4x4 = pcs->pic_disallow_4x4;

    // me_distortion/variance generated for 64x64 blocks only
    if (scs->super_block_size == 64) {
        if (rtc_tune && !pcs->ppcs->sc_class1) {
            set_depth_removal_level_controls_rtc(pcs, ctx);
        } else
            set_depth_removal_level_controls(pcs,
                                             ctx,
                                             ctx->high_freq_present
                                                 ? MAX(1,
                                                       (int)((int)pcs->pic_depth_removal_level -
                                                             (int)ctx->detect_high_freq_ctrls.depth_removal_shift))
                                                 : pcs->pic_depth_removal_level);
    }
    if (rtc_tune) {
        if (enc_mode <= ENC_M8)
            set_lpd1_ctrls(ctx,
                           ctx->high_freq_present ? MIN(pcs->pic_lpd1_lvl, ctx->detect_high_freq_ctrls.max_pic_lpd1_lvl)
                                                  : pcs->pic_lpd1_lvl);
        else {
            int lpd1_lvl = pcs->pic_lpd1_lvl;
            if (pcs->slice_type != I_SLICE) {
                int me_8x8 = pcs->ppcs->me_8x8_cost_variance[ctx->sb_index];
                int th     = enc_mode <= ENC_M11 ? 3 * ctx->qp_index : 3000;

                if (me_8x8 < th)
                    lpd1_lvl += 2;
            }
            // Checking against 6, as 6 is the max level for lpd1_lvl. If max level for lpd1_lvl changes, the check should be updated
            lpd1_lvl = MAX(0, MIN(lpd1_lvl, 6));
            set_lpd1_ctrls(ctx, lpd1_lvl);
        }
    } else
        set_lpd1_ctrls(ctx,
                       ctx->high_freq_present ? MIN(pcs->pic_lpd1_lvl, ctx->detect_high_freq_ctrls.max_pic_lpd1_lvl)
                                              : pcs->pic_lpd1_lvl);

    if (!rtc_tune)
        ctx->pd1_lvl_refinement = 0;
    else if ((pcs->ppcs->sc_class1 && enc_mode <= ENC_M10) || (!pcs->ppcs->sc_class1 && enc_mode <= ENC_M10))
        ctx->pd1_lvl_refinement = 1;
    else
        ctx->pd1_lvl_refinement = 2;
    // 1st call to avoid using invalid settings at the construction of the block(s) queue of PD1 when regular PD0 is not called
    svt_aom_set_nsq_ctrls(ctx, pcs->nsq_level, NULL, NULL, NULL);
}
static void set_depth_early_exit_ctrls(ModeDecisionContext *ctx, uint8_t early_exit_level) {
    DepthEarlyExitCtrls *ctrls = &ctx->depth_early_exit_ctrls;

    switch (early_exit_level) {
    case 0:
        ctrls->split_cost_th = 0;
        ctrls->early_exit_th = 0;
        break;
    case 1:
        ctrls->split_cost_th = 50;
        ctrls->early_exit_th = 0;
        break;
    case 2:
        ctrls->split_cost_th = 50;
        ctrls->early_exit_th = 900;
        break;

    default: assert(0); break;
    }
}
// Set signals used for light-pd0 path; only PD0 should call this function
// assumes NSQ OFF, no 4x4, no chroma, no TXT/TXS/RDOQ/SSSE, SB_64x64
void svt_aom_sig_deriv_enc_dec_light_pd0(SequenceControlSet *scs, PictureControlSet *pcs, ModeDecisionContext *ctx) {
    const Pd0Level           pd0_level = ctx->lpd0_ctrls.pd0_level;
    PictureParentControlSet *ppcs      = pcs->ppcs;
    const uint8_t            is_ref    = ppcs->is_ref;
    const uint8_t            is_islice = pcs->slice_type == I_SLICE;
    const bool               rtc_tune  = (scs->static_config.pred_structure == SVT_AV1_PRED_LOW_DELAY_B) ? true : false;

    ctx->md_disallow_nsq                 = 1;
    ctx->svt_aom_inject_inter_candidates = 1;

    // Use coeff rate and slit flag rate only (i.e. no fast rate)
    ctx->shut_fast_rate = TRUE;

    uint8_t depth_early_exit_lvl = 1;
    // When only the predicted depth is used, use safe early exit THs
    if (rtc_tune && pd0_level == VERY_LIGHT_PD0)
        depth_early_exit_lvl = 0;
    else if (pd0_level <= LPD0_LVL_0 || ctx->pic_pred_depth_only)
        depth_early_exit_lvl = 1;
    else
        depth_early_exit_lvl = 2;
    set_depth_early_exit_ctrls(ctx, depth_early_exit_lvl);

    uint8_t intra_level = 0;
    if (pcs->slice_type == I_SLICE || ppcs->transition_present == 1)
        intra_level = 1;
    else if (pd0_level <= LPD0_LVL_1) {
        uint64_t use_intra_pd0_th = compute_intra_pd0_th(scs, pcs, ctx);
        uint32_t fast_lambda      = ctx->hbd_md ? ctx->fast_lambda_md[EB_10_BIT_MD] : ctx->fast_lambda_md[EB_8_BIT_MD];
        uint64_t cost_64x64       = RDCOST(fast_lambda, 0, ppcs->me_64x64_distortion[ctx->sb_index]);

        intra_level = (cost_64x64 < use_intra_pd0_th) ? 0 : 1;
    } else
        intra_level = 0;
    set_intra_ctrls(pcs, ctx, intra_level);
    if (pd0_level == VERY_LIGHT_PD0) {
        // Modulate the inter-depth bias based on the QP and the temporal complexity of the SB
        // towards more split for low QPs or/and complex SBs,
        // and less split for high QPs or/and easy SBs (to compensate for the absence of the coeff rate)
        const uint32_t init_bias_tab[INPUT_SIZE_COUNT] = {800, 850, 850, 850, 900, 950, 1000};
        const uint32_t init_bias_value                 = (rtc_tune ? init_bias_tab[pcs->ppcs->input_resolution] : 950);
        ctx->inter_depth_bias = init_bias_value + pcs->ppcs->frm_hdr.quantization_params.base_q_idx;
        if (ppcs->me_8x8_cost_variance[ctx->sb_index] > 1000)
            ctx->inter_depth_bias = ctx->inter_depth_bias - 150;
        else if (ppcs->me_8x8_cost_variance[ctx->sb_index] > 500)
            ctx->inter_depth_bias = ctx->inter_depth_bias - 50;
        else if (ppcs->me_8x8_cost_variance[ctx->sb_index] > 250)
            ctx->inter_depth_bias = ctx->inter_depth_bias + 50;
        else
            ctx->inter_depth_bias = ctx->inter_depth_bias + 150;
    } else {
        ctx->inter_depth_bias = 0;
    }

    ctx->d2_parent_bias = 1000;
    set_mds0_controls(ctx, 4);
    if (pd0_level == VERY_LIGHT_PD0)
        return;
    svt_aom_set_chroma_controls(ctx, 0 /*chroma off*/);

    // Using PF in LPD0 may cause some VQ issues
    set_pf_controls(ctx, 1);

    uint8_t subres_level;
    // LPD0 was designed assuming 4x4 blocks were disallowed. Since LPD0 is now used in some presets where 4x4 is on
    // check that subres is not used when 4x4 blocks are enabled.
    if (pd0_level <= LPD0_LVL_1 || !ctx->disallow_4x4) {
        subres_level = 0;
    } else {
        subres_level = 0;

        // Use b64_geom here because LPD0 is only enabled when 64x64 SB size is used
        B64Geom *b64_geom = &ppcs->b64_geom[ctx->sb_index];

        // The controls checks the deviation between: (1) the pred-to-src SAD of even rows and (2) the pred-to-src SAD of odd rows for each 64x64 to decide whether to use subres or not
        // then applies the result to the 64x64 block and to all children, therefore if incomplete 64x64 then shut subres
        if (b64_geom->is_complete_b64) {
            // Use ME distortion and variance detector to enable subres
            uint64_t use_subres_th = compute_subres_th(scs, pcs, ctx);
            uint32_t fast_lambda   = ctx->hbd_md ? ctx->fast_lambda_md[EB_10_BIT_MD] : ctx->fast_lambda_md[EB_8_BIT_MD];
            uint64_t cost_64x64    = RDCOST(fast_lambda, 0, ppcs->me_64x64_distortion[ctx->sb_index]);

            if (pd0_level <= LPD0_LVL_2) {
                if (is_islice || ppcs->transition_present == 1)
                    subres_level = 1;
                else
                    subres_level = (cost_64x64 < use_subres_th) ? 1 : 0;
            } else if (pd0_level <= LPD0_LVL_3) {
                if (is_islice || ppcs->transition_present == 1)
                    subres_level = 1;
                else if (is_ref)
                    subres_level = (cost_64x64 < use_subres_th) ? 1 : 0;
                else
                    subres_level = 2;
            } else {
                if (is_ref)
                    subres_level = (ctx->depth_removal_ctrls.enabled &&
                                    (ctx->depth_removal_ctrls.disallow_below_16x16 ||
                                     ctx->depth_removal_ctrls.disallow_below_32x32 ||
                                     ctx->depth_removal_ctrls.disallow_below_64x64))
                        ? 2
                        : 1;
                else
                    subres_level = 2;
            }
        } else {
            if (pd0_level <= LPD0_LVL_2)
                subres_level = 0;
            else
                subres_level = is_ref ? 0 : 2;
        }
    }
    set_subres_controls(ctx, subres_level);

    if (pd0_level <= LPD0_LVL_2)
        set_rate_est_ctrls(ctx, 2);
    else if (pd0_level <= LPD0_LVL_3)
        set_rate_est_ctrls(ctx, 4);
    else
        set_rate_est_ctrls(ctx, 0);

    // set at pic-level b/c feature depends on some pic-level initializations
    ctx->approx_inter_rate = 1;
}

void svt_aom_sig_deriv_enc_dec_light_pd1(PictureControlSet *pcs, ModeDecisionContext *ctx) {
    Pd1Level                 lpd1_level       = ctx->lpd1_ctrls.pd1_level;
    PictureParentControlSet *ppcs             = pcs->ppcs;
    const uint8_t            is_ref           = ppcs->is_ref;
    const EbInputResolution  input_resolution = ppcs->input_resolution;
    const uint8_t            is_islice        = pcs->slice_type == I_SLICE;
    const SliceType          slice_type       = pcs->slice_type;
    // Get ref info, used to set some feature levels
    const uint32_t picture_qp           = pcs->picture_qp;
    uint32_t       me_8x8_cost_variance = (uint32_t)~0;
    uint32_t       me_64x64_distortion  = (uint32_t)~0;
    uint8_t        l0_was_skip = 0, l1_was_skip = 0;
    uint8_t        l0_was_64x64_mvp = 0, l1_was_64x64_mvp = 0;
    const bool     rtc_tune  = (pcs->scs->static_config.pred_structure == SVT_AV1_PRED_LOW_DELAY_B) ? true : false;
    const uint8_t  sc_class1 = ppcs->sc_class1;
    const EncMode  enc_mode  = pcs->enc_mode;
    const Bool     disable_rdoq_rtc = enc_mode <= ENC_M12 ? 0 : rtc_tune && sc_class1 ? 1 : 0;
    // the frame size of reference pics are different if enable reference scaling.
    // sb info can not be reused because super blocks are mismatched, so we set
    // the reference pic unavailable to avoid using wrong info
    const Bool is_ref_l0_avail = svt_aom_is_ref_same_size(pcs, REF_LIST_0, 0);
    const Bool is_ref_l1_avail = svt_aom_is_ref_same_size(pcs, REF_LIST_1, 0);

    // REF info only available if frame is not an I_SLICE
    if (!is_islice && is_ref_l0_avail) {
        me_8x8_cost_variance          = ppcs->me_8x8_cost_variance[ctx->sb_index];
        me_64x64_distortion           = ppcs->me_64x64_distortion[ctx->sb_index];
        EbReferenceObject *ref_obj_l0 = (EbReferenceObject *)pcs->ref_pic_ptr_array[REF_LIST_0][0]->object_ptr;
        l0_was_skip = ref_obj_l0->sb_skip[ctx->sb_index], l1_was_skip = 1;
        l0_was_64x64_mvp = ref_obj_l0->sb_64x64_mvp[ctx->sb_index], l1_was_64x64_mvp = 1;
        if (slice_type == B_SLICE && is_ref_l1_avail) {
            EbReferenceObject *ref_obj_l1 = (EbReferenceObject *)pcs->ref_pic_ptr_array[REF_LIST_1][0]->object_ptr;
            l1_was_skip                   = ref_obj_l1->sb_skip[ctx->sb_index];
            l1_was_64x64_mvp              = ref_obj_l1->sb_64x64_mvp[ctx->sb_index];
        }
    }
    uint8_t ref_skip_perc = pcs->ref_skip_percentage;

    // Set candidate reduction levels
    uint8_t cand_reduction_level = 0;
    if (is_islice) {
        cand_reduction_level = 0;
    } else if (rtc_tune) {
        if (lpd1_level <= LPD1_LVL_0)
            cand_reduction_level = 3;
        else if (lpd1_level <= LPD1_LVL_2)
            cand_reduction_level = 4;
        else if (lpd1_level <= LPD1_LVL_3)
            cand_reduction_level = 5;
        else
            cand_reduction_level = 6;
    } else {
        if (lpd1_level <= LPD1_LVL_0)
            cand_reduction_level = 4;
        else if (lpd1_level <= LPD1_LVL_2)
            cand_reduction_level = 5;
        else if (lpd1_level <= LPD1_LVL_3)
            cand_reduction_level = 6;
        else
            cand_reduction_level = 7;
    }
    if (ppcs->scs->rc_stat_gen_pass_mode)
        cand_reduction_level = 7;
    set_cand_reduction_ctrls(pcs,
                             ctx,
                             cand_reduction_level,
                             picture_qp,
                             me_8x8_cost_variance,
                             me_64x64_distortion,
                             l0_was_skip,
                             l1_was_skip,
                             ref_skip_perc);
    // Rdoq is disabled in this case for light pd1 regardless lpd1 level
    if (disable_rdoq_rtc) {
        ctx->rdoq_level = 0;
    } else {
        if (lpd1_level <= LPD1_LVL_0)
            ctx->rdoq_level = 4;
        else if (lpd1_level <= LPD1_LVL_4)
            ctx->rdoq_level = 5;
        else
            ctx->rdoq_level = 0;
    }
    set_rdoq_controls(ctx, ctx->rdoq_level);
    if (lpd1_level <= LPD1_LVL_0)
        ctx->md_subpel_me_level = input_resolution <= INPUT_SIZE_1080p_RANGE ? (is_ref ? 10 : 12) : 16;
    else if (lpd1_level <= LPD1_LVL_1)
        ctx->md_subpel_me_level = input_resolution <= INPUT_SIZE_1080p_RANGE ? (is_ref ? 13 : 14) : 16;
    else if (lpd1_level <= LPD1_LVL_2) {
        ctx->md_subpel_me_level = input_resolution <= INPUT_SIZE_1080p_RANGE ? (is_ref ? 14 : 15) : 16;
    } else {
        ctx->md_subpel_me_level = input_resolution <= INPUT_SIZE_1080p_RANGE ? (is_ref ? 14 : 15) : 16;
        if (((l0_was_skip && l1_was_skip && ref_skip_perc > 50) || (l0_was_64x64_mvp && l1_was_64x64_mvp)) &&
            me_8x8_cost_variance < (200 * picture_qp) && me_64x64_distortion < (200 * picture_qp))
            ctx->md_subpel_me_level = 0;
    }
    md_subpel_me_controls(ctx, ctx->md_subpel_me_level, rtc_tune);

    uint8_t mds0_level = 0;
    if (lpd1_level <= LPD1_LVL_4)
        mds0_level = 4;
    else {
        mds0_level = is_ref ? 4 : 0;
        if (((l0_was_skip && l1_was_skip && ref_skip_perc > 40) ||
             (me_8x8_cost_variance < (250 * picture_qp) && me_64x64_distortion < (250 * picture_qp))))
            mds0_level = 0;
    }

    set_mds0_controls(ctx, mds0_level);

    uint8_t lpd1_tx_level = 0;
    if (lpd1_level <= LPD1_LVL_2)
        lpd1_tx_level = 3;
    else {
        lpd1_tx_level = 4;
        if (!rtc_tune &&
            (((l0_was_skip && l1_was_skip && ref_skip_perc > 35) && me_8x8_cost_variance < (800 * picture_qp) &&
              me_64x64_distortion < (800 * picture_qp)) ||
             (me_8x8_cost_variance < (100 * picture_qp) && me_64x64_distortion < (100 * picture_qp))))
            lpd1_tx_level = 6;
    }
    set_lpd1_tx_ctrls(pcs, ctx, lpd1_tx_level);

    /* In modes below M13, only use level 1-3 for chroma detector, as more aggressive levels will cause
    blurring artifacts in certain clips.

    Do not test this signal in M12 and below during preset tuning.  This signal should be kept as an enc_mode check
    instead of and LPD1_LEVEL check to ensure that M12 and below do not use it.
    */
    if ((rtc_tune && pcs->enc_mode >= ENC_M12) || (!rtc_tune && pcs->enc_mode >= ENC_M13)) {
        if (lpd1_level <= LPD1_LVL_4)
            ctx->lpd1_tx_ctrls.chroma_detector_level = 4;
        else
            ctx->lpd1_tx_ctrls.chroma_detector_level = 0;
    }

    /* In modes below M13, only skip non-NEAREST_NEAREST TX b/c skipping all inter TX will cause blocking artifacts
    in certain clips.  This signal is separated from the general lpd1_tx_ctrls (above) to avoid
    accidentally turning this on for modes below M13.

    Do not test this signal in M12 and below during preset tuning.  This signal should be kept as an enc_mode check
    instead of and LPD1_LEVEL check to ensure that M12 and below do not use it.
    */
    if (rtc_tune) {
        if ((!pcs->ppcs->sc_class1 && pcs->enc_mode <= ENC_M11) || (pcs->ppcs->sc_class1 && pcs->enc_mode <= ENC_M10))
            ctx->lpd1_skip_inter_tx_level = 0;
        else {
            assert(pcs->enc_mode >= ENC_M11 && "Only enable this feature for M13+ in RA or M11+ for low delay");
            if (lpd1_level <= LPD1_LVL_2) {
                ctx->lpd1_skip_inter_tx_level = 0;
            } else if (lpd1_level <= LPD1_LVL_4) {
                ctx->lpd1_skip_inter_tx_level = 1;
            } else {
                ctx->lpd1_skip_inter_tx_level = 1;
                if (((l0_was_skip && l1_was_skip && ref_skip_perc > 35) && me_8x8_cost_variance < (800 * picture_qp) &&
                     me_64x64_distortion < (800 * picture_qp)) ||
                    (me_8x8_cost_variance < (100 * picture_qp) && me_64x64_distortion < (100 * picture_qp))) {
                    ctx->lpd1_skip_inter_tx_level = 2;
                }
            }
        }
    } else {
        if (pcs->enc_mode <= ENC_M12)
            ctx->lpd1_skip_inter_tx_level = 0;
        else {
            assert(pcs->enc_mode >= ENC_M13 && "Only enable this feature for M13+ in RA or M11+ for low delay");
            if (lpd1_level <= LPD1_LVL_2) {
                ctx->lpd1_skip_inter_tx_level = 0;
            } else if (lpd1_level <= LPD1_LVL_4) {
                ctx->lpd1_skip_inter_tx_level = 1;
                if (((l0_was_skip && l1_was_skip && ref_skip_perc > 35) && me_8x8_cost_variance < (800 * picture_qp) &&
                     me_64x64_distortion < (800 * picture_qp)) ||
                    (me_8x8_cost_variance < (100 * picture_qp) && me_64x64_distortion < (100 * picture_qp))) {
                    ctx->lpd1_skip_inter_tx_level = 2;
                }
            } else {
                ctx->lpd1_skip_inter_tx_level = is_ref ? 1 : 2;
                if (((l0_was_skip && l1_was_skip && ref_skip_perc > 35) && me_8x8_cost_variance < (800 * picture_qp) &&
                     me_64x64_distortion < (800 * picture_qp)) ||
                    (me_8x8_cost_variance < (100 * picture_qp) && me_64x64_distortion < (100 * picture_qp))) {
                    ctx->lpd1_skip_inter_tx_level = 2;
                }
            }
        }
    }
    // 0: Feature off
    // Lower the threshold, the more aggressive the feature is
    ctx->lpd1_bypass_tx_th_div = 0;
    if (pcs->rtc_tune && !pcs->ppcs->sc_class1) {
        ctx->lpd1_bypass_tx_th_div = enc_mode <= ENC_M8 ? 0 : enc_mode <= ENC_M10 ? 8 : 6;
    }
    uint8_t rate_est_level = 0;
    if (lpd1_level <= LPD1_LVL_0)
        rate_est_level = 4;
    else
        rate_est_level = 0;
    set_rate_est_ctrls(ctx, rate_est_level);

    // If want to turn off approximating inter rate, must ensure that the approximation is also disabled
    // at the pic level (pcs->approx_inter_rate)
    ctx->approx_inter_rate = 1;

    uint8_t pf_level = 1;
    if (lpd1_level <= LPD1_LVL_4)
        pf_level = 1;
    else
        pf_level = 2;
    set_pf_controls(ctx, pf_level);

    uint8_t intra_level = 0;
    if (lpd1_level <= LPD1_LVL_2)
        intra_level = 6;
    else
        intra_level = 7;
    set_intra_ctrls(pcs, ctx, intra_level);
    ctx->d2_parent_bias = 995;
    /* Set signals that have assumed values in the light-PD1 path (but need to be initialized as they may be checked) */

    // Use coeff rate and slit flag rate only (i.e. no fast rate)
    ctx->shut_fast_rate   = FALSE;
    ctx->uv_ctrls.enabled = 1;
    ctx->uv_ctrls.uv_mode = CHROMA_MODE_1;
    set_cfl_ctrls(ctx, 0);
    ctx->md_disallow_nsq                            = 1;
    ctx->new_nearest_injection                      = 1;
    ctx->svt_aom_inject_inter_candidates            = 1;
    ctx->blk_skip_decision                          = TRUE;
    ctx->rate_est_ctrls.update_skip_ctx_dc_sign_ctx = 0;
    ctx->rate_est_ctrls.update_skip_coeff_ctx       = 0;
    ctx->subres_ctrls.odd_to_even_deviation_th      = 0;
}
void svt_aom_sig_deriv_enc_dec(SequenceControlSet *scs, PictureControlSet *pcs, ModeDecisionContext *ctx) {
    EncMode                  enc_mode             = pcs->enc_mode;
    uint8_t                  pd_pass              = ctx->pd_pass;
    PictureParentControlSet *ppcs                 = pcs->ppcs;
    const uint8_t            is_ref               = ppcs->is_ref;
    const uint8_t            is_base              = ppcs->temporal_layer_index == 0;
    const EbInputResolution  input_resolution     = ppcs->input_resolution;
    const uint8_t            is_islice            = pcs->slice_type == I_SLICE;
    const uint32_t           hierarchical_levels  = ppcs->hierarchical_levels;
    const uint32_t           picture_qp           = pcs->picture_qp;
    uint32_t                 me_8x8_cost_variance = (uint32_t)~0;
    uint32_t                 me_64x64_distortion  = (uint32_t)~0;
    uint8_t                  l0_was_skip = 0, l1_was_skip = 0;
    uint8_t                  ref_skip_perc = pcs->ref_skip_percentage;

    // 2nd call as set_nsq_ctrls() has a PD_PASS check
    svt_aom_set_nsq_ctrls(ctx, pcs->nsq_level, NULL, NULL, NULL);

    const bool rtc_tune = (scs->static_config.pred_structure == SVT_AV1_PRED_LOW_DELAY_B) ? true : false;

    set_cand_reduction_ctrls(pcs,
                             ctx,
                             pd_pass == PD_PASS_0 ? 0 : pcs->cand_reduction_level,
                             picture_qp,
                             me_8x8_cost_variance,
                             me_64x64_distortion,
                             l0_was_skip,
                             l1_was_skip,
                             ref_skip_perc);

    uint8_t txt_level = pcs->txt_level;
    if (ctx->high_freq_present)
        txt_level = txt_level == 0 ? ctx->detect_high_freq_ctrls.max_pd1_txt_lvl
                                   : MIN(txt_level, ctx->detect_high_freq_ctrls.max_pd1_txt_lvl);
    svt_aom_set_txt_controls(ctx, pd_pass == PD_PASS_0 ? 0 : txt_level);
    set_tx_shortcut_ctrls(pcs, ctx, pd_pass == PD_PASS_0 ? 0 : ctx->high_freq_present ? 0 : pcs->tx_shortcut_level);

    set_interpolation_search_level_ctrls(ctx, pd_pass == PD_PASS_0 ? 0 : pcs->interpolation_search_level);

    svt_aom_set_chroma_controls(ctx, pd_pass == PD_PASS_0 ? 0 : pcs->chroma_level);

    set_cfl_ctrls(ctx, pd_pass == PD_PASS_0 ? 0 : pcs->cfl_level);
    if (pd_pass == PD_PASS_0)
        ctx->md_disallow_nsq = enc_mode <= ENC_MR ? !ctx->nsq_ctrls.enabled : 1;
    else {
        // Update nsq settings based on the sb_class
        ctx->md_disallow_nsq = !ctx->nsq_ctrls.enabled;
    }

    if (pd_pass == PD_PASS_0)
        ctx->global_mv_injection = 0;
    else
        ctx->global_mv_injection = ppcs->gm_ctrls.enabled;
    if (pd_pass == PD_PASS_0)
        ctx->new_nearest_injection = 0;
    else
        ctx->new_nearest_injection = 1;
    ctx->new_nearest_near_comb_injection = pd_pass == PD_PASS_0 ? 0 : pcs->new_nearest_near_comb_injection;

    //set Warped-Motion controls from Picture level.
    svt_aom_set_wm_controls(ctx, pd_pass == PD_PASS_0 ? 0 : pcs->wm_level);

    ctx->unipred3x3_injection            = pd_pass == PD_PASS_0 ? 0 : pcs->unipred3x3_injection;
    ctx->bipred3x3_injection             = pd_pass == PD_PASS_0 ? 0 : pcs->bipred3x3_injection;
    ctx->svt_aom_inject_inter_candidates = 1;
    ctx->inter_compound_mode             = pd_pass == PD_PASS_0 ? 0 : pcs->inter_compound_mode;
    svt_aom_set_dist_based_ref_pruning_controls(ctx, pd_pass == PD_PASS_0 ? 0 : pcs->dist_based_ref_pruning);
    set_spatial_sse_full_loop_level(ctx, pd_pass == PD_PASS_0 ? 0 : pcs->spatial_sse_full_loop_level);
    if (ctx->uv_ctrls.uv_mode <= CHROMA_MODE_1)
        ctx->blk_skip_decision = TRUE;
    else
        ctx->blk_skip_decision = FALSE;
    if (pd_pass == PD_PASS_0) {
        if (enc_mode <= ENC_M3)
            ctx->rdoq_level = 1;
        else
            ctx->rdoq_level = 0;
    } else if (rtc_tune) {
        if (ppcs->sc_class1) {
            if (enc_mode <= ENC_M11)
                ctx->rdoq_level = 1;
            else
                ctx->rdoq_level = 0;
        } else
            ctx->rdoq_level = 1;
    } else {
        if (enc_mode <= ENC_M12)
            ctx->rdoq_level = 1;
        else
            ctx->rdoq_level = 5;
    }
    if (scs->low_latency_kf && is_islice)
        ctx->rdoq_level = 0;
    set_rdoq_controls(ctx, ctx->rdoq_level);

    // Derive redundant block
    if (pd_pass == PD_PASS_0 || ctx->md_disallow_nsq)
        ctx->redundant_blk = FALSE;
    else
        ctx->redundant_blk = TRUE;
    set_parent_sq_coeff_area_based_cycles_reduction_ctrls(
        ctx, input_resolution, pd_pass == PD_PASS_0 ? 0 : ctx->nsq_ctrls.psq_cplx_lvl);
    set_sq_txs_ctrls(ctx, ctx->nsq_ctrls.psq_txs_lvl);
    set_sq_pred_ctrls(ctx, ctx->nsq_ctrls.psq_pred_lvl);
    uint8_t depth_early_exit_lvl = 0;
    if (pd_pass == PD_PASS_0)
        depth_early_exit_lvl = 1;
    else if (enc_mode <= ENC_M7)
        depth_early_exit_lvl = 1;
    else
        depth_early_exit_lvl = 2;
    set_depth_early_exit_ctrls(ctx, depth_early_exit_lvl);
    // Set pic_obmc_level @ MD
    if (pd_pass == PD_PASS_0)
        ctx->md_pic_obmc_level = 0;
    else
        ctx->md_pic_obmc_level = pcs->ppcs->pic_obmc_level;

    svt_aom_set_obmc_controls(ctx, ctx->md_pic_obmc_level);
    svt_aom_set_inter_intra_ctrls(ctx, pd_pass == PD_PASS_0 ? 0 : pcs->inter_intra_level);
    set_txs_controls(ctx, pd_pass == PD_PASS_0 ? 0 : pcs->txs_level);
    // Set md_filter_intra_mode @ MD
    // md_filter_intra_level specifies whether filter intra would be active
    // for a given prediction candidate in mode decision.

    // md_filter_intra_level | Settings
    // 0                      | OFF
    // 1                      | ON
    if (pd_pass == PD_PASS_0)
        ctx->md_filter_intra_level = 0;
    else
        ctx->md_filter_intra_level = pcs->pic_filter_intra_level;
    // Set md_allow_intrabc @ MD
    if (pd_pass == PD_PASS_0)
        ctx->md_allow_intrabc = 0;
    else
        ctx->md_allow_intrabc = pcs->ppcs->frm_hdr.allow_intrabc;

    // Set md_palette_level @ MD
    if (pd_pass == PD_PASS_0)
        ctx->md_palette_level = 0;
    else
        ctx->md_palette_level = pcs->ppcs->palette_level;

    uint8_t pf_level = 1;
    set_pf_controls(ctx, pf_level);
    md_sq_motion_search_controls(ctx, pd_pass == PD_PASS_0 ? 0 : pcs->md_sq_mv_search_level);

    md_nsq_motion_search_controls(ctx, pd_pass == PD_PASS_0 ? 0 : pcs->md_nsq_mv_search_level);

    svt_aom_md_pme_search_controls(ctx, pd_pass == PD_PASS_0 ? 0 : pcs->md_pme_level);
    if (pd_pass == PD_PASS_0)
        ctx->md_subpel_me_level = enc_mode <= ENC_M5 ? 4 : 0;
    else if (enc_mode <= ENC_MR)
        ctx->md_subpel_me_level = 1;
    else if (enc_mode <= ENC_M2)
        ctx->md_subpel_me_level = input_resolution <= INPUT_SIZE_480p_RANGE ? 1 : 2;
    else if (enc_mode <= ENC_M4)
        ctx->md_subpel_me_level = 2;
    else if (enc_mode <= ENC_M7)
        ctx->md_subpel_me_level = 3;
    else if (enc_mode <= ENC_M10)
        ctx->md_subpel_me_level = is_base ? 3 : (is_ref ? 5 : 8);
    else
        ctx->md_subpel_me_level = is_ref ? 5 : 8;
    ctx->md_subpel_me_level += hierarchical_levels <= 2 ? 2 : 0;
    ctx->md_subpel_me_level = MIN(16, ctx->md_subpel_me_level);
    md_subpel_me_controls(ctx, ctx->md_subpel_me_level, rtc_tune);
    if (pd_pass == PD_PASS_0)
        ctx->md_subpel_pme_level = enc_mode <= ENC_MR ? 3 : 0;
    else if (enc_mode <= ENC_MR)
        ctx->md_subpel_pme_level = 1;
    else
        ctx->md_subpel_pme_level = 2;

    md_subpel_pme_controls(ctx, ctx->md_subpel_pme_level);
    uint8_t rate_est_level = 0;
    if (pd_pass == PD_PASS_0) {
        if (enc_mode <= ENC_M2)
            rate_est_level = 1;
        else
            rate_est_level = 2;
    } else {
        if (enc_mode <= ENC_M12)
            rate_est_level = 1;
        else if ((enc_mode <= ENC_M12) || (rtc_tune && enc_mode <= ENC_M13))
            rate_est_level = is_islice ? 3 : 4;
        else
            rate_est_level = is_islice ? 3 : 0;
    }
    set_rate_est_ctrls(ctx, rate_est_level);

    // set at pic-level b/c feature depends on some pic-level initializations
    ctx->approx_inter_rate = pcs->approx_inter_rate;
    // Use coeff rate and slit flag rate only (i.e. no fast rate)
    if (pd_pass == PD_PASS_0)
        ctx->shut_fast_rate = TRUE;
    else
        ctx->shut_fast_rate = FALSE;

    // intra_level must be greater than 0 for I_SLICE
    uint8_t intra_level = 0;
    if (pd_pass == PD_PASS_0) {
        intra_level = 7;
    } else if (rtc_tune) {
        if (enc_mode <= ENC_M9)
            intra_level = (is_islice || ppcs->transition_present == 1) ? 1 : 6;
        else
            intra_level = (is_islice || ppcs->transition_present == 1) ? 4 : 6;
    } else if (enc_mode <= ENC_MR)
        intra_level = 1;
    else if (enc_mode <= ENC_M0)
        intra_level = is_base ? 1 : 2;
    else if (enc_mode <= ENC_M2)
        intra_level = is_base ? 1 : 3;
    else if (enc_mode <= ENC_M3)
        intra_level = is_base ? 1 : 4;
    else if (enc_mode <= ENC_M5)
        intra_level = is_base ? 1 : 5;
    else if (enc_mode <= ENC_M9)
        intra_level = (is_islice || ppcs->transition_present == 1) ? 1 : is_base ? 2 : 6;
    else
        intra_level = (is_islice || ppcs->transition_present == 1) ? 4 : 6;

    if (pcs->scs->low_latency_kf && is_islice)
        intra_level = 6;
    set_intra_ctrls(pcs, ctx, intra_level);

    set_mds0_controls(ctx, pd_pass == PD_PASS_0 ? 2 : pcs->mds0_level);

    set_subres_controls(ctx, 0);
    ctx->inter_depth_bias = 0;
    if (pd_pass == PD_PASS_0)
        ctx->d2_parent_bias = 1000;
    else
        ctx->d2_parent_bias = 995;
    uint8_t skip_sub_depth_lvl;
    if (pd_pass == PD_PASS_0 || pcs->ppcs->sc_class1)
        skip_sub_depth_lvl = 0;
    else if (enc_mode <= ENC_M0)
        skip_sub_depth_lvl = 1;
    else if (enc_mode <= ENC_M5)
        skip_sub_depth_lvl = 2;
    else
        skip_sub_depth_lvl = 3;

    set_skip_sub_depth_ctrls(&ctx->skip_sub_depth_ctrls, skip_sub_depth_lvl);
    ctx->tune_ssim_level = SSIM_LVL_0;
}
/*
* return the 4x4 level
Used by svt_aom_sig_deriv_enc_dec and memory allocation
*/
bool svt_aom_get_disallow_4x4(EncMode enc_mode, SliceType slice_type) {
    (void)slice_type;
    if (enc_mode <= ENC_M5)
        return false;
    else
        return true;
}
// Get the nsq_level used for each preset (to be passed to setting function: svt_aom_set_nsq_ctrls())
uint8_t svt_aom_get_nsq_level(EncMode enc_mode, uint8_t is_islice, uint8_t is_base, InputCoeffLvl coeff_lvl) {
    uint8_t nsq_level;
    //set the nsq_level
    if (enc_mode <= ENC_M0)
        nsq_level = is_islice ? 2 : 3;
    else if (enc_mode <= ENC_M1)
        nsq_level = is_base ? 5 : 6;
    else if (enc_mode <= ENC_M4)
        nsq_level = is_base ? 7 : 8;
    else if (enc_mode <= ENC_M6) {
        if (coeff_lvl == LOW_LVL)
            nsq_level = is_base ? 17 : 20;
        else if (coeff_lvl == HIGH_LVL)
            nsq_level = is_base ? 19 : 22;
        else // regular
            nsq_level = is_base ? 18 : 21;
    } else if (enc_mode <= ENC_M7) {
        if (coeff_lvl == LOW_LVL)
            nsq_level = 29;
        else if (coeff_lvl == HIGH_LVL)
            nsq_level = 31;
        else // regular
            nsq_level = 30;
    } else
        nsq_level = 0;

    return nsq_level;
}
/*
* return by-pass encdec
*/
uint8_t svt_aom_get_bypass_encdec(EncMode enc_mode, uint8_t hbd_md, uint8_t encoder_bit_depth) {
    UNUSED(hbd_md);
    uint8_t bypass_encdec = 1;
    if (encoder_bit_depth == EB_EIGHT_BIT) {
        // 8bit settings
        if (enc_mode <= ENC_M5)
            bypass_encdec = 0;
        else
            bypass_encdec = 1;
    } else {
        // 10bit settings
        // If the preset is changed for this feature, the warning in svt_av1_verify_settings in
        // EbEncSetting.c must also be changed to reflect the change in preset
        // This feature causes a mismatch between the files encoded for 10bit when bypass_encdec
        // is 1 if recon is enabled/disabled or stat-report is enabled/disabled
        if (enc_mode <= ENC_M9)
            bypass_encdec = 0;
        else
            bypass_encdec = 1;
    }
    return bypass_encdec;
}

// use this function to set the disallow_below_16x16 level in MD. ME 8x8 blocks are controlled by svt_aom_get_enable_me_8x8()
static uint8_t svt_aom_get_disallow_below_16x16_picture_level(EncMode enc_mode, EbInputResolution resolution,
                                                              Bool is_islice, Bool sc_class1, Bool is_ref,
                                                              bool rtc_tune) {
    uint8_t disallow_below_16x16 = 0;
    if (sc_class1)
        disallow_below_16x16 = 0;
    else if ((enc_mode <= ENC_M8) || (rtc_tune && (enc_mode <= ENC_M9)))
        disallow_below_16x16 = 0;
    else if (enc_mode <= ENC_M13) {
        if (resolution <= INPUT_SIZE_1080p_RANGE)
            disallow_below_16x16 = is_ref ? 0 : 1;
        else
            disallow_below_16x16 = is_islice ? 0 : 1;
    } else
        disallow_below_16x16 = is_islice ? 0 : 1;

    return disallow_below_16x16;
}

static void set_cdf_controls(PictureControlSet *pcs, uint8_t update_cdf_level) {
    CdfControls *ctrl = &pcs->cdf_ctrl;
    switch (update_cdf_level) {
    case 0:
        ctrl->update_mv   = 0;
        ctrl->update_se   = 0;
        ctrl->update_coef = 0;
        break;
    case 1:
        ctrl->update_mv   = 1;
        ctrl->update_se   = 1;
        ctrl->update_coef = 1;
        break;
    case 2:
        ctrl->update_mv   = 0;
        ctrl->update_se   = 1;
        ctrl->update_coef = 1;
        break;
    case 3:
        ctrl->update_mv   = 0;
        ctrl->update_se   = 1;
        ctrl->update_coef = 0;
        break;
    default: assert(0); break;
    }

    ctrl->update_mv = pcs->slice_type == I_SLICE ? 0 : ctrl->update_mv;
    ctrl->enabled   = ctrl->update_coef | ctrl->update_mv | ctrl->update_se;
}
/******************************************************
* Derive Mode Decision Config Settings for OQ
Input   : encoder mode and tune
Output  : EncDec Kernel signal(s)
******************************************************/
static EbErrorType rtime_alloc_ec_ctx_array(PictureControlSet *pcs, uint16_t all_sb) {
    EB_MALLOC_ARRAY(pcs->ec_ctx_array, all_sb);
    return EB_ErrorNone;
}

uint8_t svt_aom_get_update_cdf_level(EncMode enc_mode, SliceType is_islice, uint8_t is_base) {
    uint8_t update_cdf_level = 0;
    if (enc_mode <= ENC_M2)
        update_cdf_level = 1;
    else if (enc_mode <= ENC_M5)
        update_cdf_level = is_base ? 1 : 3;
    else if (enc_mode <= ENC_M9)
        update_cdf_level = is_islice ? 1 : 0;
    else
        update_cdf_level = 0;

    return update_cdf_level;
}

uint8_t svt_aom_get_chroma_level(EncMode enc_mode) {
    uint8_t chroma_level = 0;
    if (enc_mode <= ENC_MR)
        chroma_level = 2;
    else if (enc_mode <= ENC_M12)
        chroma_level = 4;
    else
        chroma_level = 5;

    return chroma_level;
}

/*
set lpd0_level
*/
static void set_pic_lpd0_lvl(PictureControlSet *pcs, EncMode enc_mode) {
    PictureParentControlSet *ppcs = pcs->ppcs;

    const uint8_t is_base            = ppcs->temporal_layer_index == 0;
    const uint8_t is_islice          = pcs->slice_type == I_SLICE;
    const Bool    transition_present = (ppcs->transition_present == 1);
    const bool    rtc_tune  = (ppcs->scs->static_config.pred_structure == SVT_AV1_PRED_LOW_DELAY_B) ? true : false;
    InputCoeffLvl coeff_lvl = pcs->coeff_lvl;
    const EbInputResolution input_resolution = ppcs->input_resolution;

    if (rtc_tune) {
        if (enc_mode <= ENC_M9) {
            if (coeff_lvl == LOW_LVL) {
                pcs->pic_lpd0_lvl = 2;
            } else if (coeff_lvl == HIGH_LVL) {
                pcs->pic_lpd0_lvl = (is_base || transition_present) ? 3 : 5;
            } else { // Regular
                pcs->pic_lpd0_lvl = 3;
            }
        } else if (enc_mode <= ENC_M10 && !pcs->ppcs->sc_class1) {
            if (coeff_lvl == LOW_LVL) {
                pcs->pic_lpd0_lvl = 3;
            } else if (coeff_lvl == HIGH_LVL) {
                pcs->pic_lpd0_lvl = (is_base || transition_present) ? 6 : 7;
            } else { // Regular
                pcs->pic_lpd0_lvl = (is_base || transition_present) ? 3 : 5;
            }
        } else if (enc_mode <= ENC_M11) {
            if (coeff_lvl == LOW_LVL) {
                pcs->pic_lpd0_lvl = 3;
            } else if (coeff_lvl == HIGH_LVL) {
                pcs->pic_lpd0_lvl = (is_base || transition_present) ? 6 : 8;
            } else { // Regular
                pcs->pic_lpd0_lvl = (is_base || transition_present) ? 5 : 7;
            }
        } else if (enc_mode <= ENC_M12) {
            if (coeff_lvl == LOW_LVL) {
                pcs->pic_lpd0_lvl = 5;
            } else if (coeff_lvl == HIGH_LVL) {
                pcs->pic_lpd0_lvl = (is_base || transition_present) ? 6 : 8;
            } else { // Regular
                pcs->pic_lpd0_lvl = (is_base || transition_present) ? 5 : 7;
            }
        } else
            pcs->pic_lpd0_lvl = is_islice ? 5 : 8;
    } else {
        if (enc_mode <= ENC_M3)
            pcs->pic_lpd0_lvl = 0;
        else if (enc_mode <= ENC_M8)
            pcs->pic_lpd0_lvl = 1;
        else if (enc_mode <= ENC_M9) {
            if (coeff_lvl == LOW_LVL) {
                pcs->pic_lpd0_lvl = 2;
            } else if (coeff_lvl == HIGH_LVL) {
                pcs->pic_lpd0_lvl = (is_base || transition_present) ? 3 : 5;
            } else { // Regular
                pcs->pic_lpd0_lvl = 3;
            }
        } else if (enc_mode <= ENC_M10) {
            if (input_resolution <= INPUT_SIZE_360p_RANGE)
                pcs->pic_lpd0_lvl = 3;
            else if (input_resolution <= INPUT_SIZE_480p_RANGE)
                pcs->pic_lpd0_lvl = (is_base || transition_present) ? 3 : 5;
            else {
                if (pcs->coeff_lvl == HIGH_LVL)
                    pcs->pic_lpd0_lvl = (is_base || transition_present) ? 5 : 7;
                else
                    pcs->pic_lpd0_lvl = (is_base || transition_present) ? 3 : 5;
            }
        } else if (enc_mode <= ENC_M12) {
            if (coeff_lvl == LOW_LVL) {
                pcs->pic_lpd0_lvl = (is_base || transition_present) ? 3 : 5;
            } else if (coeff_lvl == HIGH_LVL) {
                pcs->pic_lpd0_lvl = (is_base || transition_present) ? 6 : 8;
            } else { // Regular
                pcs->pic_lpd0_lvl = (is_base || transition_present) ? 6 : 7;
            }
        } else {
            if (coeff_lvl == LOW_LVL) {
                pcs->pic_lpd0_lvl = (is_base || transition_present) ? 6 : 7;
            } else if (coeff_lvl == HIGH_LVL) {
                pcs->pic_lpd0_lvl = 8;
            } else { // regular
                pcs->pic_lpd0_lvl = (is_base || transition_present) ? 6 : 8;
            }
        }
    }
}
uint8_t get_inter_compound_level(EncMode enc_mode) {
    uint8_t comp_level;
    if (enc_mode <= ENC_M0)
        comp_level = 3;
    else if (enc_mode <= ENC_M2)
        comp_level = 4;
    else
        comp_level = 0;

    return comp_level;
}
uint8_t get_filter_intra_level(EncMode enc_mode) {
    uint8_t filter_intra_level;
    if (enc_mode <= ENC_M3)
        filter_intra_level = 1;
    else
        filter_intra_level = 0;

    return filter_intra_level;
}
uint8_t svt_aom_get_inter_intra_level(EncMode enc_mode, uint8_t is_base, uint8_t transition_present) {
    uint8_t inter_intra_level = 0;
    if (enc_mode <= ENC_M1)
        inter_intra_level = 2;
    else if (enc_mode <= ENC_M2)
        inter_intra_level = (transition_present || is_base) ? 2 : 0;
    else if (enc_mode <= ENC_M11)
        inter_intra_level = transition_present ? 2 : 0;
    else
        inter_intra_level = 0;

    return inter_intra_level;
}

uint8_t svt_aom_get_obmc_level(EncMode enc_mode, uint8_t fast_decode, EbInputResolution input_resolution) {
    uint8_t obmc_level = 0;
    if (fast_decode == 0 || input_resolution <= INPUT_SIZE_360p_RANGE) {
        if (enc_mode <= ENC_M1)
            obmc_level = 1;
        else if (enc_mode <= ENC_M2)
            obmc_level = 2;
        else if (enc_mode <= ENC_M6)
            obmc_level = 3;
        else if (enc_mode <= ENC_M7)
            obmc_level = 4;
        else
            obmc_level = 0;
    } else {
        if (enc_mode <= ENC_M1)
            obmc_level = 1;
        else if (enc_mode <= ENC_M7)
            obmc_level = 3;
        else
            obmc_level = 0;
    }
    return obmc_level;
}
void svt_aom_sig_deriv_mode_decision_config(SequenceControlSet *scs, PictureControlSet *pcs) {
    PictureParentControlSet *ppcs                = pcs->ppcs;
    EncMode                  enc_mode            = pcs->enc_mode;
    const uint8_t            is_ref              = ppcs->is_ref;
    const uint8_t            is_base             = ppcs->temporal_layer_index == 0;
    const uint8_t            is_layer1           = ppcs->temporal_layer_index == 1;
    const EbInputResolution  input_resolution    = ppcs->input_resolution;
    const uint8_t            is_islice           = pcs->slice_type == I_SLICE;
    const SliceType          slice_type          = pcs->slice_type;
    const Bool               fast_decode         = scs->static_config.fast_decode;
    const uint32_t           hierarchical_levels = ppcs->hierarchical_levels;
    const Bool               transition_present  = (ppcs->transition_present == 1);
    const bool               rtc_tune = (scs->static_config.pred_structure == SVT_AV1_PRED_LOW_DELAY_B) ? true : false;
    //MFMV
    if (is_islice || scs->mfmv_enabled == 0 || pcs->ppcs->frm_hdr.error_resilient_mode) {
        ppcs->frm_hdr.use_ref_frame_mvs = 0;
    } else {
        if (fast_decode == 0) {
            if ((!rtc_tune && enc_mode <= ENC_M9) || (rtc_tune && (enc_mode <= ENC_M8)))
                ppcs->frm_hdr.use_ref_frame_mvs = 1;
            else {
                uint64_t avg_me_dist = 0;
                for (uint16_t b64_idx = 0; b64_idx < ppcs->b64_total_count; b64_idx++) {
                    avg_me_dist += ppcs->me_64x64_distortion[b64_idx];
                }
                avg_me_dist /= ppcs->b64_total_count;
                avg_me_dist /= pcs->picture_qp;

                ppcs->frm_hdr.use_ref_frame_mvs = avg_me_dist < 200 || input_resolution <= INPUT_SIZE_360p_RANGE ? 1
                                                                                                                 : 0;
            }
        } else {
            if (enc_mode <= ENC_M10) {
                uint64_t avg_me_dist = 0;
                for (uint16_t b64_idx = 0; b64_idx < ppcs->b64_total_count; b64_idx++) {
                    avg_me_dist += ppcs->me_64x64_distortion[b64_idx];
                }
                avg_me_dist /= ppcs->b64_total_count;
                avg_me_dist /= pcs->picture_qp;

                ppcs->frm_hdr.use_ref_frame_mvs = avg_me_dist < 50 || input_resolution <= INPUT_SIZE_360p_RANGE ? 1 : 0;
            } else {
                ppcs->frm_hdr.use_ref_frame_mvs = input_resolution <= INPUT_SIZE_360p_RANGE ? 1 : 0;
            }
        }
    }

    uint8_t update_cdf_level = svt_aom_get_update_cdf_level(enc_mode, is_islice, is_base);
    //set the conrols uisng the required level
    set_cdf_controls(pcs, update_cdf_level);

    if (pcs->cdf_ctrl.enabled) {
        const uint16_t picture_sb_w = ppcs->picture_sb_width;
        const uint16_t picture_sb_h = ppcs->picture_sb_height;
        const uint16_t all_sb       = picture_sb_w * picture_sb_h;
        rtime_alloc_ec_ctx_array(pcs, all_sb);
    }
    //Filter Intra Mode : 0: OFF  1: ON
    // pic_filter_intra_level specifies whether filter intra would be active
    // for a given picture.
    pcs->pic_filter_intra_level = get_filter_intra_level(enc_mode);
    if (fast_decode == 0 || input_resolution <= INPUT_SIZE_360p_RANGE) {
        if (pcs->enc_mode <= ENC_M6)
            pcs->ppcs->use_accurate_part_ctx = TRUE;
        else
            pcs->ppcs->use_accurate_part_ctx = FALSE;
    } else {
        pcs->ppcs->use_accurate_part_ctx = FALSE;
    }
    FrameHeader *frm_hdr             = &ppcs->frm_hdr;
    frm_hdr->allow_high_precision_mv = frm_hdr->quantization_params.base_q_idx < HIGH_PRECISION_MV_QTHRESH &&
            (scs->input_resolution <= INPUT_SIZE_480p_RANGE)
        ? 1
        : 0;
    // Set Warped Motion level and enabled flag
    pcs->wm_level = 0;
    if (frm_hdr->frame_type == KEY_FRAME || frm_hdr->frame_type == INTRA_ONLY_FRAME || frm_hdr->error_resilient_mode ||
        pcs->ppcs->frame_superres_enabled || pcs->ppcs->frame_resize_enabled) {
        pcs->wm_level = 0;
    } else {
        if (fast_decode == 0 || input_resolution <= INPUT_SIZE_360p_RANGE) {
            if (enc_mode <= ENC_M0) {
                pcs->wm_level = 1;
            } else if (enc_mode <= ENC_M3) {
                pcs->wm_level = 2;

            } else if (enc_mode <= ENC_M7) {
                if (hierarchical_levels <= 3)
                    pcs->wm_level = is_base ? 1 : 0;
                else
                    pcs->wm_level = (is_base || is_layer1) ? 2 : 0;
            } else if (enc_mode <= ENC_M8) {
                pcs->wm_level = (is_base || is_layer1) ? 2 : 0;
            } else if (enc_mode <= ENC_M13) {
                if (input_resolution <= INPUT_SIZE_720p_RANGE)
                    pcs->wm_level = is_base ? 3 : 0;
                else
                    pcs->wm_level = is_base ? 4 : 0;
            } else {
                pcs->wm_level = is_base ? 4 : 0;
            }
        } else {
            if (enc_mode <= ENC_M7) {
                pcs->wm_level = is_base ? 1 : 0;
            } else if (enc_mode <= ENC_M9) {
                if (input_resolution <= INPUT_SIZE_720p_RANGE)
                    pcs->wm_level = is_base ? 1 : 0;
                else
                    pcs->wm_level = is_base ? 4 : 0;
            } else {
                pcs->wm_level = is_base ? 4 : 0;
            }
        }
    }
    if (hierarchical_levels <= 2) {
        pcs->wm_level = enc_mode <= ENC_M7 ? pcs->wm_level : 0;
    }
    Bool enable_wm = pcs->wm_level ? 1 : 0;
    if (scs->enable_warped_motion != DEFAULT)
        enable_wm = (Bool)scs->enable_warped_motion;
    // Note: local warp should be disabled when super-res or resize is ON
    // according to the AV1 spec 5.11.27
    frm_hdr->allow_warped_motion = enable_wm &&
        !(frm_hdr->frame_type == KEY_FRAME || frm_hdr->frame_type == INTRA_ONLY_FRAME) &&
        !frm_hdr->error_resilient_mode && !pcs->ppcs->frame_superres_enabled &&
        scs->static_config.resize_mode == RESIZE_NONE;

    frm_hdr->is_motion_mode_switchable = frm_hdr->allow_warped_motion;
    ppcs->pic_obmc_level               = svt_aom_get_obmc_level(enc_mode, fast_decode, input_resolution);
    // Switchable Motion Mode
    frm_hdr->is_motion_mode_switchable = frm_hdr->is_motion_mode_switchable || ppcs->pic_obmc_level;

    ppcs->bypass_cost_table_gen = 0;
    if (rtc_tune) {
        if (enc_mode <= ENC_M10)
            pcs->approx_inter_rate = 0;
        else
            pcs->approx_inter_rate = 1;
    } else {
        if (enc_mode <= ENC_M11)
            pcs->approx_inter_rate = 0;
        else
            pcs->approx_inter_rate = 1;
    }
    if (is_islice || transition_present)
        pcs->skip_intra = 0;
    if (rtc_tune) {
        if (enc_mode <= ENC_M9)
            pcs->skip_intra = 0;
        else
            pcs->skip_intra = (is_ref || pcs->ref_intra_percentage > 50) ? 0 : 1;
    } else {
        if (enc_mode <= ENC_M8)
            pcs->skip_intra = 0;
        else
            pcs->skip_intra = (is_ref || pcs->ref_intra_percentage > 50) ? 0 : 1;
    }

    // Set the level for the candidate(s) reduction feature
    pcs->cand_reduction_level = 0;
    if (is_islice)
        pcs->cand_reduction_level = 0;
    else if (enc_mode <= ENC_MR)
        pcs->cand_reduction_level = 0;
    else if (enc_mode <= ENC_M3)
        pcs->cand_reduction_level = is_base ? 0 : 1;
    else if (enc_mode <= ENC_M7) {
        if (pcs->coeff_lvl == LOW_LVL)
            pcs->cand_reduction_level = 1;
        else
            pcs->cand_reduction_level = 2;
    } else if (rtc_tune) {
        if (enc_mode <= ENC_M9) {
            if (pcs->coeff_lvl == LOW_LVL)
                pcs->cand_reduction_level = 1;
            else
                pcs->cand_reduction_level = 2;
        } else
            pcs->cand_reduction_level = 2;
    } else {
        if (enc_mode <= ENC_M12) {
            if (pcs->coeff_lvl == LOW_LVL)
                pcs->cand_reduction_level = 1;
            else
                pcs->cand_reduction_level = 3;
        } else
            pcs->cand_reduction_level = 3;
    }
    if (scs->rc_stat_gen_pass_mode)
        pcs->cand_reduction_level = 7;
    // Set the level for the txt search
    pcs->txt_level = 0;
    if (rtc_tune) {
        if ((!pcs->ppcs->sc_class1 && (enc_mode <= ENC_M11)) || (pcs->ppcs->sc_class1 && (enc_mode <= ENC_M13))) {
            if (pcs->coeff_lvl == LOW_LVL) {
                pcs->txt_level = is_base ? 4 : 5;
            } else if (pcs->coeff_lvl == HIGH_LVL) {
                pcs->txt_level = is_base ? 8 : 10;
            } else { // regular
                pcs->txt_level = 7;
            }
        } else {
            if (pcs->coeff_lvl == LOW_LVL) {
                pcs->txt_level = is_base ? 8 : 10;
            } else if (pcs->coeff_lvl == HIGH_LVL) {
                pcs->txt_level = is_base ? 8 : 10;
                if (pcs->ref_intra_percentage < 85 && !pcs->ppcs->sc_class1) {
                    pcs->txt_level = 0;
                }
            } else { // regular
                pcs->txt_level = is_base ? 8 : 10;
                if (pcs->ref_intra_percentage < 85 && !is_base && !pcs->ppcs->sc_class1) {
                    pcs->txt_level = 0;
                }
            }
        }
    } else {
        if (enc_mode <= ENC_M1)
            pcs->txt_level = 2;
        else if (enc_mode <= ENC_M7) {
            if (pcs->coeff_lvl == LOW_LVL) {
                pcs->txt_level = 4;
            } else if (pcs->coeff_lvl == HIGH_LVL) {
                pcs->txt_level = 7;
            } else { // regular
                pcs->txt_level = is_base ? 4 : 5;
            }
        } else if (enc_mode <= ENC_M12) {
            if (pcs->coeff_lvl == LOW_LVL) {
                pcs->txt_level = is_base ? 4 : 5;
            } else if (pcs->coeff_lvl == HIGH_LVL) {
                pcs->txt_level = is_base ? 8 : 10;
            } else { // regular
                pcs->txt_level = 7;
            }
        } else {
            pcs->txt_level = is_base ? 8 : 10;
            if (pcs->ref_intra_percentage < 85 && !pcs->ppcs->sc_class1) {
                pcs->txt_level = 0;
            }
        }
    }
    // Set the level for the txt shortcut feature
    // Any tx_shortcut_level having the chroma detector off in REF frames should be reserved for M13+
    pcs->tx_shortcut_level = 0;
    if (rtc_tune) {
        if (enc_mode <= ENC_M4)
            pcs->tx_shortcut_level = 0;
        else if ((enc_mode <= ENC_M7 && !pcs->ppcs->sc_class1) || (enc_mode <= ENC_M8 && pcs->ppcs->sc_class1))
            pcs->tx_shortcut_level = is_base ? 0 : 1;
        else if (enc_mode <= ENC_M10)
            pcs->tx_shortcut_level = is_islice ? 0 : 1;
        else
            pcs->tx_shortcut_level = is_islice ? 0 : 4;
    } else {
        if (enc_mode <= ENC_M5)
            pcs->tx_shortcut_level = 0;

        else if (enc_mode <= ENC_M13)
            pcs->tx_shortcut_level = is_base ? 0 : 1;
        else
            pcs->tx_shortcut_level = is_islice ? 0 : 4;
    }
    // Set the level the interpolation search
    pcs->interpolation_search_level = 0;
    if (enc_mode <= ENC_M5)
        pcs->interpolation_search_level = 4;
    else {
        pcs->interpolation_search_level = 4;
        if (!is_base) {
            const uint8_t th[INPUT_SIZE_COUNT] = {100, 100, 85, 50, 30, 30, 30};
            const uint8_t skip_area            = pcs->ref_skip_percentage;
            if (skip_area > th[input_resolution])
                pcs->interpolation_search_level = 0;
        }
    }
    frm_hdr->interpolation_filter = SWITCHABLE;

    pcs->chroma_level = svt_aom_get_chroma_level(enc_mode);
    // Set the level for cfl
    pcs->cfl_level = 0;
    if (pcs->ppcs->sc_class1) {
        if (enc_mode <= ENC_M6)
            pcs->cfl_level = 1;
        else
            pcs->cfl_level = is_base ? 2 : 0;
    } else if (enc_mode <= ENC_M3)
        pcs->cfl_level = 1;
    else if (rtc_tune) {
        if (enc_mode <= ENC_M8) {
            pcs->cfl_level = is_base ? 2 : 0;
        } else if (enc_mode <= ENC_M11) {
            if (hierarchical_levels <= 3)
                pcs->cfl_level = is_islice ? 2 : 0;
            else
                pcs->cfl_level = is_base ? 2 : 0;
        } else if (enc_mode <= ENC_M13)
            pcs->cfl_level = is_islice ? 2 : 0;
        else
            pcs->cfl_level = 0;
    } else {
        if (enc_mode <= ENC_M9) {
            pcs->cfl_level = is_base ? 2 : 0;
        } else if (enc_mode <= ENC_M12)
            pcs->cfl_level = is_islice ? 2 : 0;
        else
            pcs->cfl_level = 0;
    }
    if (scs->low_latency_kf && is_islice)
        pcs->cfl_level = 0;
    // Set the level for new/nearest/near injection
    if (scs->new_nearest_comb_inject == DEFAULT)
        if (enc_mode <= ENC_MR)
            pcs->new_nearest_near_comb_injection = 1;
        else
            pcs->new_nearest_near_comb_injection = 0;
    else
        pcs->new_nearest_near_comb_injection = scs->new_nearest_comb_inject;

    // Set the level for unipred3x3 injection
    if (enc_mode <= ENC_MR)
        pcs->unipred3x3_injection = 1;
    else
        pcs->unipred3x3_injection = 0;

    // Set the level for bipred3x3 injection
    if (scs->bipred_3x3_inject == DEFAULT) {
        if (enc_mode <= ENC_MR)
            pcs->bipred3x3_injection = 1;
        else if (enc_mode <= ENC_M2)
            pcs->bipred3x3_injection = 2;
        else
            pcs->bipred3x3_injection = 0;
    } else {
        pcs->bipred3x3_injection = scs->bipred_3x3_inject;
    }

    // Set the level for inter-inter compound
    pcs->inter_compound_mode = get_inter_compound_level(enc_mode);

    // Set the level for the distance-based red pruning
    if (pcs->ppcs->ref_list0_count_try > 1 || pcs->ppcs->ref_list1_count_try > 1) {
        if (enc_mode <= ENC_M4)
            pcs->dist_based_ref_pruning = is_base ? 2 : 5;
        else if (!rtc_tune && enc_mode <= ENC_M11)
            pcs->dist_based_ref_pruning = is_base ? 2 : 6;
        else {
            if (pcs->coeff_lvl == LOW_LVL) {
                pcs->dist_based_ref_pruning = is_base ? 2 : 6;
            } else {
                pcs->dist_based_ref_pruning = is_base ? 3 : 6;
            }
        }
    } else {
        pcs->dist_based_ref_pruning = 0;
    }

    // Set the level the spatial sse @ full-loop
    pcs->spatial_sse_full_loop_level = 0;
    if (scs->spatial_sse_full_loop_level == DEFAULT)
        if (pcs->ppcs->sc_class1)
            pcs->spatial_sse_full_loop_level = 1;
        else if ((enc_mode <= ENC_M11) || (rtc_tune && enc_mode <= ENC_M13))
            pcs->spatial_sse_full_loop_level = 1;
        else
            pcs->spatial_sse_full_loop_level = 0;
    else
        pcs->spatial_sse_full_loop_level = scs->spatial_sse_full_loop_level;
    //set the nsq_level
    pcs->nsq_level = svt_aom_get_nsq_level(enc_mode, is_islice, is_base, pcs->coeff_lvl);

    // Set the level for inter-intra level
    if (!is_islice && scs->seq_header.enable_interintra_compound) {
        pcs->inter_intra_level = svt_aom_get_inter_intra_level(enc_mode, is_base, transition_present);
    } else
        pcs->inter_intra_level = 0;
    if (enc_mode <= ENC_M2)
        pcs->txs_level = 2;
    else if (enc_mode <= ENC_M8)
        pcs->txs_level = is_base ? 2 : 0;
    else if (rtc_tune) {
        if (enc_mode <= ENC_M11)
            pcs->txs_level = is_islice ? 3 : 0;
        else
            pcs->txs_level = is_islice ? 5 : 0;
    } else {
        if (enc_mode <= ENC_M9)
            pcs->txs_level = is_islice ? 3 : 0;
        else if (enc_mode <= ENC_M11) {
            if (hierarchical_levels <= 3)
                pcs->txs_level = is_islice ? 5 : 0;
            else
                pcs->txs_level = is_islice ? 3 : 0;
        } else
            pcs->txs_level = is_islice ? 5 : 0;
    }
    if (pcs->scs->low_latency_kf && pcs->slice_type == I_SLICE)
        pcs->txs_level = 5;
    // Set tx_mode for the frame header
    frm_hdr->tx_mode = (pcs->txs_level) ? TX_MODE_SELECT : TX_MODE_LARGEST;
    // Set the level for nic
    pcs->nic_level = svt_aom_get_nic_level(enc_mode, is_base, hierarchical_levels, rtc_tune);
    // Set the level for SQ me-search
    if (enc_mode <= ENC_MR)
        pcs->md_sq_mv_search_level = 1;
    else
        pcs->md_sq_mv_search_level = 0;

    // Set the level for NSQ me-search
    pcs->md_nsq_mv_search_level = 4;
    // Set the level for PME search
    if (enc_mode <= ENC_M0)
        pcs->md_pme_level = 2;
    else if (enc_mode <= ENC_M5) {
        if (hierarchical_levels <= 3)
            pcs->md_pme_level = 7;
        else
            pcs->md_pme_level = 3;
    } else if (enc_mode <= ENC_M7) {
        if (hierarchical_levels <= 3)
            pcs->md_pme_level = 7;
        else
            pcs->md_pme_level = 5;
    } else if (enc_mode <= ENC_M9) {
        if (hierarchical_levels <= 4)
            pcs->md_pme_level = 6;
        else
            pcs->md_pme_level = 5;
    } else
        pcs->md_pme_level = 7;
    // Set the level for mds0
    pcs->mds0_level = 0;
    if (rtc_tune) {
        if (enc_mode <= ENC_M11)
            pcs->mds0_level = 2;
        else
            pcs->mds0_level = is_islice ? 2 : 4;
    } else {
        if (enc_mode <= ENC_M9)
            pcs->mds0_level = 2;
        else if (enc_mode <= ENC_M11) {
            if (hierarchical_levels <= 3)
                pcs->mds0_level = is_islice ? 2 : 4;
            else
                pcs->mds0_level = 2;
        } else
            pcs->mds0_level = is_islice ? 2 : 4;
    }
    /*
       disallow_4x4
    */
    pcs->pic_disallow_4x4 = svt_aom_get_disallow_4x4(enc_mode, slice_type);
    /*
       Bypassing EncDec
    */
    // TODO: Bypassing EncDec doesn't work if NSQ is enabled for 10bit content (causes r2r).
    // TODO: This signal can only be modified per picture right now, not per SB.  Per SB requires
    // neighbour array updates at EncDec for all SBs, that are currently skipped if EncDec is bypassed.
    // TODO: Bypassing EncDec doesn't work if pcs->cdf_ctrl.update_coef is enabled for non-ISLICE frames (causes r2r)
    if ((scs->static_config.encoder_bit_depth == EB_EIGHT_BIT || !pcs->nsq_level) &&
        (!pcs->cdf_ctrl.update_coef || is_islice) && !ppcs->frm_hdr.segmentation_params.segmentation_enabled) {
        pcs->pic_bypass_encdec = svt_aom_get_bypass_encdec(
            enc_mode, ppcs->hbd_md, scs->static_config.encoder_bit_depth);
    } else
        pcs->pic_bypass_encdec = 0;

    /*
        set lpd0_level
    */
    // for the low delay enhance base layer frames, lower the enc_mode to improve the quality
    set_pic_lpd0_lvl(pcs, (pcs->ppcs->ld_enhanced_base_frame) ? enc_mode - 1 : enc_mode);

    if (pcs->ppcs->sc_class1 || scs->static_config.pass == ENC_MIDDLE_PASS)
        pcs->pic_skip_pd0 = 0;
    else if (enc_mode <= ENC_M13)
        pcs->pic_skip_pd0 = 0;
    else
        pcs->pic_skip_pd0 = is_base ? 0 : 1;
    pcs->pic_disallow_below_16x16 = svt_aom_get_disallow_below_16x16_picture_level(
        enc_mode, input_resolution, is_islice, ppcs->sc_class1, is_ref, rtc_tune);

    if (scs->super_block_size == 64) {
        if (is_islice || transition_present) {
            pcs->pic_depth_removal_level = 0;
        } else {
            // Set depth_removal_level_controls
            if (pcs->ppcs->sc_class1) {
                if (enc_mode <= ENC_M7 || (!rtc_tune && enc_mode <= ENC_M8)) {
                    pcs->pic_depth_removal_level = 0;
                } else if (enc_mode <= ENC_M10) {
                    pcs->pic_depth_removal_level = is_base ? 0 : 6;
                } else if (enc_mode <= ENC_M12) {
                    pcs->pic_depth_removal_level = is_base ? 4 : 6;
                } else {
                    pcs->pic_depth_removal_level = is_base ? 5 : 14;
                }
            } else if (fast_decode == 0) {
                if (enc_mode <= ENC_M6)
                    pcs->pic_depth_removal_level = 0;
                else if ((rtc_tune && (enc_mode <= ENC_M9))) {
                    if (pcs->coeff_lvl == LOW_LVL) {
                        if (input_resolution <= INPUT_SIZE_480p_RANGE)
                            pcs->pic_depth_removal_level = 1;
                        else
                            pcs->pic_depth_removal_level = 2;
                    } else {
                        if (input_resolution <= INPUT_SIZE_480p_RANGE)
                            pcs->pic_depth_removal_level = 1;
                        else
                            pcs->pic_depth_removal_level = is_base ? 2 : 6;
                    }
                } else if (enc_mode <= ENC_M10) {
                    if (pcs->coeff_lvl == LOW_LVL) {
                        if (input_resolution <= INPUT_SIZE_480p_RANGE)
                            pcs->pic_depth_removal_level = 1;
                        else
                            pcs->pic_depth_removal_level = 2;
                    } else {
                        if (input_resolution <= INPUT_SIZE_360p_RANGE)
                            pcs->pic_depth_removal_level = is_base ? 2 : 3;
                        else if (input_resolution <= INPUT_SIZE_480p_RANGE)
                            pcs->pic_depth_removal_level = is_base ? 2 : 5;
                        else
                            pcs->pic_depth_removal_level = is_base ? 2 : 6;
                    }
                } else if (enc_mode <= ENC_M11) {
                    if (input_resolution <= INPUT_SIZE_360p_RANGE)
                        pcs->pic_depth_removal_level = is_base ? 2 : 3;
                    else if (input_resolution <= INPUT_SIZE_480p_RANGE)
                        pcs->pic_depth_removal_level = is_base ? 2 : 5;
                    else if (input_resolution <= INPUT_SIZE_720p_RANGE)
                        pcs->pic_depth_removal_level = is_base ? 2 : 6;
                    else if (input_resolution <= INPUT_SIZE_1080p_RANGE)
                        pcs->pic_depth_removal_level = is_base ? 3 : 8;
                    else
                        pcs->pic_depth_removal_level = is_base ? 9 : 14;
                } else {
                    if (input_resolution <= INPUT_SIZE_360p_RANGE)
                        pcs->pic_depth_removal_level = is_base ? 2 : 3;
                    else if (input_resolution <= INPUT_SIZE_480p_RANGE)
                        pcs->pic_depth_removal_level = is_base ? 9 : 11;
                    else
                        pcs->pic_depth_removal_level = is_base ? 9 : 14;
                }
            } else {
                if (enc_mode <= ENC_M2)
                    pcs->pic_depth_removal_level = 0;
                else if (enc_mode <= ENC_M5) {
                    if (input_resolution <= INPUT_SIZE_480p_RANGE)
                        pcs->pic_depth_removal_level = 1;
                    else
                        pcs->pic_depth_removal_level = 2;
                } else if (enc_mode <= ENC_M6) {
                    if (input_resolution <= INPUT_SIZE_1080p_RANGE)
                        pcs->pic_depth_removal_level = 1;
                    else
                        pcs->pic_depth_removal_level = 2;
                } else if (enc_mode <= ENC_M9) {
                    if (input_resolution <= INPUT_SIZE_360p_RANGE)
                        pcs->pic_depth_removal_level = is_base ? 2 : 3;
                    else if (input_resolution <= INPUT_SIZE_480p_RANGE)
                        pcs->pic_depth_removal_level = is_base ? 2 : 5;
                    else
                        pcs->pic_depth_removal_level = is_base ? 2 : 6;
                } else if (enc_mode <= ENC_M11) {
                    if (input_resolution <= INPUT_SIZE_360p_RANGE)
                        pcs->pic_depth_removal_level = is_base ? 2 : 4;
                    else if (input_resolution <= INPUT_SIZE_480p_RANGE)
                        pcs->pic_depth_removal_level = is_base ? 2 : 5;
                    else if (input_resolution <= INPUT_SIZE_720p_RANGE)
                        pcs->pic_depth_removal_level = is_base ? 2 : 6;
                    else if (input_resolution <= INPUT_SIZE_1080p_RANGE)
                        pcs->pic_depth_removal_level = is_base ? 3 : 8;
                    else
                        pcs->pic_depth_removal_level = is_base ? 9 : 14;
                } else {
                    if (input_resolution <= INPUT_SIZE_360p_RANGE)
                        pcs->pic_depth_removal_level = 7;
                    else if (input_resolution <= INPUT_SIZE_480p_RANGE)
                        pcs->pic_depth_removal_level = is_base ? 9 : 11;
                    else
                        pcs->pic_depth_removal_level = is_base ? 9 : 14;
                }
            }
        }
    }
    pcs->pic_depth_removal_level_rtc = enc_mode <= ENC_M8 || (enc_mode <= ENC_M10 && is_ref) ? 0 : 1;
    if (ppcs->sc_class1) {
        if (enc_mode <= ENC_M6)
            pcs->pic_block_based_depth_refinement_level = 0;
        else if (enc_mode <= ENC_M9)
            pcs->pic_block_based_depth_refinement_level = is_base ? 0 : hierarchical_levels == 5 ? 5 : 6;
        else if (enc_mode <= ENC_M10)
            pcs->pic_block_based_depth_refinement_level = hierarchical_levels == 5 ? (is_islice ? 1 : 5)
                                                                                   : (is_islice ? 2 : 6);
        else
            pcs->pic_block_based_depth_refinement_level = hierarchical_levels == 5 ? (is_islice ? 7 : 12)
                                                                                   : (is_islice ? 8 : 13);
    } else if (enc_mode <= ENC_M3)
        pcs->pic_block_based_depth_refinement_level = 0;
    else if (enc_mode <= ENC_M7) {
        if (pcs->coeff_lvl == LOW_LVL)
            pcs->pic_block_based_depth_refinement_level = is_base ? 1 : hierarchical_levels == 5 ? 2 : 3;
        else
            pcs->pic_block_based_depth_refinement_level = is_base ? 1 : hierarchical_levels == 5 ? 4 : 5;
    } else if (enc_mode <= ENC_M8) {
        if (rtc_tune) {
            if (pcs->coeff_lvl == LOW_LVL) {
                pcs->pic_block_based_depth_refinement_level = is_base ? 2 : 5;
            } else { // regular
                pcs->pic_block_based_depth_refinement_level = is_base ? 3 : 6;
            }
        } else {
            if (pcs->coeff_lvl == LOW_LVL) {
                pcs->pic_block_based_depth_refinement_level = hierarchical_levels == 5 ? (is_base ? 1 : 4)
                                                                                       : (is_base ? 2 : 5);
            } else { // regular
                pcs->pic_block_based_depth_refinement_level = hierarchical_levels == 5 ? (is_base ? 2 : 5)
                                                                                       : (is_base ? 3 : 6);
            }
        }
    } else {
        if (rtc_tune) {
            if (pcs->coeff_lvl == LOW_LVL) {
                pcs->pic_block_based_depth_refinement_level = is_base ? 3 : 7;
            } else { // regular
                pcs->pic_block_based_depth_refinement_level = is_base ? 8 : 11;
            }
        } else {
            if (pcs->coeff_lvl == LOW_LVL) {
                pcs->pic_block_based_depth_refinement_level = hierarchical_levels == 5 ? (is_base ? 2 : 6)
                                                                                       : (is_base ? 3 : 7);
            } else { // regular
                pcs->pic_block_based_depth_refinement_level = hierarchical_levels == 5 ? (is_base ? 7 : 10)
                                                                                       : (is_base ? 8 : 11);
            }
        }
    }
    if (pcs->ppcs->sc_class1) {
        if (enc_mode <= ENC_M7)
            pcs->pic_lpd1_lvl = 0;
        else if (enc_mode <= ENC_M9)
            pcs->pic_lpd1_lvl = is_ref ? 0 : 1;
        else if (enc_mode <= ENC_M10)
            pcs->pic_lpd1_lvl = is_ref ? 0 : 2;
        else if (enc_mode <= ENC_M11)
            pcs->pic_lpd1_lvl = is_base ? 0 : 2;
        else
            pcs->pic_lpd1_lvl = is_base ? 0 : 4;
    } else if (rtc_tune) {
        if (enc_mode <= ENC_M9) {
            if (pcs->coeff_lvl == LOW_LVL) {
                pcs->pic_lpd1_lvl = 0;
            } else if (pcs->coeff_lvl == HIGH_LVL) {
                pcs->pic_lpd1_lvl = is_base ? 0 : 2;
            } else { // Regular
                pcs->pic_lpd1_lvl = is_ref ? 0 : 1;
            }
        } else if (enc_mode <= ENC_M10) {
            if (pcs->coeff_lvl == LOW_LVL) {
                pcs->pic_lpd1_lvl = is_ref ? 0 : 1;
            } else if (pcs->coeff_lvl == HIGH_LVL) {
                pcs->pic_lpd1_lvl = is_base ? 0 : 3;
            } else { // Regular
                pcs->pic_lpd1_lvl = is_base ? 0 : 2;
            }
        } else if (enc_mode <= ENC_M12) {
            if (pcs->coeff_lvl == LOW_LVL) {
                pcs->pic_lpd1_lvl = is_ref ? 0 : 2;
            } else if (pcs->coeff_lvl == HIGH_LVL) {
                pcs->pic_lpd1_lvl = is_base ? 0 : 4;
            } else { // Regular
                pcs->pic_lpd1_lvl = is_base ? 0 : 3;
            }
        } else {
            pcs->pic_lpd1_lvl = is_base ? 0 : 4;
        }
    } else {
        if (enc_mode <= ENC_M7)
            pcs->pic_lpd1_lvl = 0;
        else if (enc_mode <= ENC_M8) {
            if (pcs->coeff_lvl == LOW_LVL) {
                pcs->pic_lpd1_lvl = 0;
            } else if (pcs->coeff_lvl == HIGH_LVL) {
                pcs->pic_lpd1_lvl = is_base ? 0 : 2;
            } else { // Regular
                pcs->pic_lpd1_lvl = is_ref ? 0 : 1;
            }
        } else if (enc_mode <= ENC_M10) {
            if (input_resolution <= INPUT_SIZE_360p_RANGE)
                pcs->pic_lpd1_lvl = is_ref ? 0 : 1;
            else if (input_resolution <= INPUT_SIZE_480p_RANGE)
                pcs->pic_lpd1_lvl = is_base ? 0 : 1;
            else {
                if (pcs->coeff_lvl == HIGH_LVL)
                    pcs->pic_lpd1_lvl = is_base ? 0 : 3;
                else
                    pcs->pic_lpd1_lvl = is_base ? 0 : 1;
            }
        } else if (enc_mode <= ENC_M12) {
            if (pcs->coeff_lvl == LOW_LVL) {
                pcs->pic_lpd1_lvl = is_base ? 0 : 2;
            } else if (pcs->coeff_lvl == HIGH_LVL) {
                pcs->pic_lpd1_lvl = is_base ? 0 : 4;
            } else { // Regular
                pcs->pic_lpd1_lvl = is_base ? 0 : 3;
            }
        } else {
            if (input_resolution <= INPUT_SIZE_1080p_RANGE && scs->static_config.encoder_bit_depth == EB_EIGHT_BIT)
                pcs->pic_lpd1_lvl = is_base ? 0 : 6;
            else
                pcs->pic_lpd1_lvl = is_base ? 0 : 5;
        }
    }
    if (is_base && !pcs->ppcs->ld_enhanced_base_frame && pcs->slice_type != I_SLICE &&
        ((!pcs->ppcs->sc_class1 && enc_mode > ENC_M10) || (pcs->ppcs->sc_class1 && enc_mode > ENC_M11)) &&
        pcs->rtc_tune)
        pcs->pic_lpd1_lvl = pcs->pic_lpd1_lvl + 1;
    // Can only use light-PD1 under the following conditions
    // There is another check before PD1 is called; pred_depth_only is not checked here, because some modes
    // may force pred_depth_only at the light-pd1 detector
    if (pcs->pic_lpd1_lvl && !(ppcs->hbd_md == 0 && pcs->pic_disallow_4x4 == TRUE && scs->super_block_size == 64)) {
        pcs->pic_lpd1_lvl = 0;
    }

    // Use the me-SAD-to-SATD deviation (of the 32x32 blocks) to detect the presence of isolated edges.
    // An SB is tagged as problematic when the deviation is higher than the normal (i.e. when me-sad and satd are not correlated)
    // For the detected SB(s), apply a better level for Depth-removal, LPD0, LPD1, and TXT of regular PD1.
    // Not applicable for I_SLICE and for SB 128x128
    if (pcs->slice_type == I_SLICE || scs->super_block_size == 128) {
        pcs->vq_ctrls.detect_high_freq_lvl = 0;
    } else {
        if (rtc_tune) {
            if (pcs->ppcs->sc_class1) {
                if (enc_mode <= ENC_M9)
                    pcs->vq_ctrls.detect_high_freq_lvl = 1;
                else
                    pcs->vq_ctrls.detect_high_freq_lvl = 0;
            } else
                pcs->vq_ctrls.detect_high_freq_lvl = 2;
        } else {
            if (pcs->ppcs->sc_class1) {
                pcs->vq_ctrls.detect_high_freq_lvl = 1;
            } else if (enc_mode <= ENC_M11) {
                pcs->vq_ctrls.detect_high_freq_lvl = 2;
            } else {
                pcs->vq_ctrls.detect_high_freq_lvl = 0;
            }
        }
    }
}
