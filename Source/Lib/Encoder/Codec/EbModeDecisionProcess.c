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

#include "EbUtility.h"
#include "EbModeDecisionProcess.h"
#include "EbLambdaRateTables.h"
#include "EbRateControlProcess.h"
#include "EncModeConfig.h"

void set_block_based_depth_refinement_controls(ModeDecisionContext *ctx, uint8_t block_based_depth_refinement_level);
static void mode_decision_context_dctor(EbPtr p) {
    ModeDecisionContext *obj = (ModeDecisionContext *)p;

    uint32_t block_max_count_sb = obj->init_max_block_cnt;
    for (int cd = 0; cd < MAX_PAL_CAND; cd++)
        if (obj->palette_cand_array[cd].color_idx_map)
            EB_FREE_ARRAY(obj->palette_cand_array[cd].color_idx_map);

    // MD palette search
    if (obj->palette_size_array_0)
        EB_FREE_ARRAY(obj->palette_size_array_0);
    for (CandClass cand_class_it = CAND_CLASS_0; cand_class_it < CAND_CLASS_TOTAL; cand_class_it++)
        EB_FREE_ARRAY(obj->cand_buff_indices[cand_class_it]);
    EB_FREE_ARRAY(obj->best_candidate_index_array);

    EB_FREE_ARRAY(obj->above_txfm_context);
    EB_FREE_ARRAY(obj->left_txfm_context);
    for (uint32_t coded_leaf_index = 0; coded_leaf_index < block_max_count_sb; ++coded_leaf_index) {
        if (obj->md_blk_arr_nsq[coded_leaf_index].coeff_tmp)
            EB_DELETE(obj->md_blk_arr_nsq[coded_leaf_index].coeff_tmp);
        if (obj->md_blk_arr_nsq[coded_leaf_index].recon_tmp)
            EB_DELETE(obj->md_blk_arr_nsq[coded_leaf_index].recon_tmp);
    }
    EB_DELETE_PTR_ARRAY(obj->cand_bf_ptr_array, obj->max_nics_uv);
    EB_FREE_ARRAY(obj->cand_bf_tx_depth_1->cand);
    EB_DELETE(obj->cand_bf_tx_depth_1);
    EB_FREE_ARRAY(obj->cand_bf_tx_depth_2->cand);
    EB_DELETE(obj->cand_bf_tx_depth_2);
    EB_FREE_ALIGNED_ARRAY(obj->cfl_temp_luma_recon16bit);
    EB_FREE_ALIGNED_ARRAY(obj->cfl_temp_luma_recon);
    EB_FREE_ARRAY(obj->fast_cand_array);
    EB_FREE_ARRAY(obj->fast_cand_ptr_array);
    EB_FREE_2D(obj->injected_mvs);
    EB_FREE_ARRAY(obj->injected_ref_types);
    EB_FREE_ARRAY(obj->fast_cost_array);
    EB_FREE_ARRAY(obj->full_cost_array);
    if (obj->md_local_blk_unit) {
        for (int i = 0; i < 3; i++) {
            EB_FREE_ARRAY(obj->md_local_blk_unit[0].neigh_left_recon_16bit[i]);
            EB_FREE_ARRAY(obj->md_local_blk_unit[0].neigh_top_recon_16bit[i]);
            EB_FREE_ARRAY(obj->md_local_blk_unit[0].neigh_left_recon[i]);
            EB_FREE_ARRAY(obj->md_local_blk_unit[0].neigh_top_recon[i]);
        }
    }
    if (obj->md_blk_arr_nsq) {
        EB_FREE_ARRAY(obj->md_blk_arr_nsq[0].av1xd);
    }
    EB_FREE_ARRAY(obj->avail_blk_flag);
    EB_FREE_ARRAY(obj->cost_avail);
    EB_FREE_ARRAY(obj->md_local_blk_unit);
    EB_FREE_ARRAY(obj->md_blk_arr_nsq);
    if (obj->rate_est_table)
        EB_FREE_ARRAY(obj->rate_est_table);
    if (obj->pred0)
        EB_FREE(obj->pred0);
    if (obj->pred1)
        EB_FREE(obj->pred1);
    if (obj->residual1)
        EB_FREE(obj->residual1);
    if (obj->diff10)
        EB_FREE(obj->diff10);

    if (obj->intrapred_buf)
        EB_FREE_2D(obj->intrapred_buf);

    if (obj->obmc_buff_0)
        EB_FREE(obj->obmc_buff_0);
    if (obj->obmc_buff_1)
        EB_FREE(obj->obmc_buff_1);
    if (obj->wsrc_buf)
        EB_FREE(obj->wsrc_buf);
    if (obj->mask_buf)
        EB_FREE(obj->mask_buf);
    EB_FREE_ARRAY(obj->mdc_sb_array);
    for (uint32_t txt_itr = 0; txt_itr < TX_TYPES; ++txt_itr) {
        EB_DELETE(obj->recon_coeff_ptr[txt_itr]);
        EB_DELETE(obj->recon_ptr[txt_itr]);
        EB_DELETE(obj->quant_coeff_ptr[txt_itr]);
    }
    EB_DELETE(obj->tx_coeffs);
    EB_DELETE(obj->scratch_prediction_ptr);
    EB_DELETE(obj->temp_residual);
    EB_DELETE(obj->temp_recon_ptr);
    EB_FREE_ARRAY(obj->full_cost_ssim_array);
}

void svt_aom_set_nics(NicScalingCtrls *scaling_ctrls, uint32_t mds1_count[CAND_CLASS_TOTAL],
                      uint32_t mds2_count[CAND_CLASS_TOTAL], uint32_t mds3_count[CAND_CLASS_TOTAL], uint8_t pic_type);

/******************************************************
 * Mode Decision Context Constructor
 ******************************************************/
EbErrorType svt_aom_mode_decision_context_ctor(ModeDecisionContext *ctx, EbColorFormat color_format, uint8_t sb_size,
                                               EncMode enc_mode, uint16_t max_block_cnt, uint32_t encoder_bit_depth,
                                               EbFifo *mode_decision_configuration_input_fifo_ptr,
                                               EbFifo *mode_decision_output_fifo_ptr, uint8_t enable_hbd_mode_decision,
                                               uint8_t cfg_palette, bool rtc_tune) {
    uint32_t buffer_index;
    uint32_t cand_index;

    ctx->init_max_block_cnt     = max_block_cnt;
    uint32_t block_max_count_sb = max_block_cnt;

    ctx->sb_size = sb_size;
    (void)color_format;

    ctx->dctor  = mode_decision_context_dctor;
    ctx->hbd_md = enable_hbd_mode_decision;

    // Input/Output System Resource Manager FIFOs
    ctx->mode_decision_configuration_input_fifo_ptr = mode_decision_configuration_input_fifo_ptr;
    ctx->mode_decision_output_fifo_ptr              = mode_decision_output_fifo_ptr;

    // Maximum number of candidates MD can support
    // determine MAX_NICS for a given preset
    // get the min scaling level (the smallest scaling level is the most conservative)
    uint8_t min_nic_scaling_level = NICS_SCALING_LEVELS - 1;
    for (uint8_t hl = 0; hl < MAX_TEMPORAL_LAYERS; hl++) {
        for (uint8_t is_base = 0; is_base < 2; is_base++) {
            uint8_t nic_level         = svt_aom_get_nic_level(enc_mode, is_base, hl, rtc_tune);
            uint8_t nic_scaling_level = svt_aom_set_nic_controls(NULL, nic_level);
            min_nic_scaling_level     = MIN(min_nic_scaling_level, nic_scaling_level);
        }
    }
    uint8_t stage1_scaling_num = MD_STAGE_NICS_SCAL_NUM[min_nic_scaling_level][MD_STAGE_1];

    // scale max_nics
    uint32_t max_nics = 0;
    {
        NicScalingCtrls scaling_ctrls;
        scaling_ctrls.stage1_scaling_num = stage1_scaling_num;
        scaling_ctrls.stage2_scaling_num = stage1_scaling_num;
        scaling_ctrls.stage3_scaling_num = stage1_scaling_num;
        uint32_t mds1_count[CAND_CLASS_TOTAL];
        uint32_t mds2_count[CAND_CLASS_TOTAL];
        uint32_t mds3_count[CAND_CLASS_TOTAL];
        for (uint8_t pic_type = 0; pic_type < NICS_PIC_TYPE; pic_type++) {
            svt_aom_set_nics(&scaling_ctrls, mds1_count, mds2_count, mds3_count, pic_type);

            uint32_t nics = 0;
            for (CandClass cidx = CAND_CLASS_0; cidx < CAND_CLASS_TOTAL; cidx++) { nics += mds1_count[cidx]; }
            max_nics = MAX(max_nics, nics);
        }
    }

    // If independent chroma search is used, need to allocate additional 84 candidate buffers
    const uint8_t ind_uv_cands = svt_aom_set_chroma_controls(NULL, svt_aom_get_chroma_level(enc_mode)) == CHROMA_MODE_0
        ? 84
        : 0;
    max_nics += CAND_CLASS_TOTAL; //need one extra temp buffer for each fast loop call
    ctx->max_nics    = max_nics;
    ctx->max_nics_uv = max_nics + ind_uv_cands;
    // Cfl scratch memory
    if (ctx->hbd_md > EB_8_BIT_MD)
        EB_MALLOC_ALIGNED(ctx->cfl_temp_luma_recon16bit, sizeof(uint16_t) * sb_size * sb_size);
    if (ctx->hbd_md != EB_10_BIT_MD)
        EB_MALLOC_ALIGNED(ctx->cfl_temp_luma_recon, sizeof(uint8_t) * sb_size * sb_size);
    uint8_t use_update_cdf = 0;
    for (uint8_t is_islice = 0; is_islice < 2; is_islice++) {
        for (uint8_t is_base = 0; is_base < 2; is_base++) {
            if (use_update_cdf)
                break;
            use_update_cdf |= svt_aom_get_update_cdf_level(enc_mode, is_islice, is_base);
        }
    }
    if (use_update_cdf)
        EB_CALLOC_ARRAY(ctx->rate_est_table, 1);
    else
        ctx->rate_est_table = NULL;
    // Allocate buffer for inter-inter compound prediction
    if (get_inter_compound_level(enc_mode)) {
        const uint8_t bits = ctx->hbd_md > EB_8_BIT_MD ? 2 : 1;
        EB_MALLOC(ctx->pred0, sb_size * sb_size * bits * sizeof(ctx->pred0[0]));
        EB_MALLOC(ctx->pred1, sb_size * sb_size * bits * sizeof(ctx->pred1[0]));
        EB_MALLOC(ctx->residual1, sb_size * sb_size * sizeof(ctx->residual1[0]));
        EB_MALLOC(ctx->diff10, sb_size * sb_size * sizeof(ctx->diff10[0]));
    }

    // Allocate buffer for inter-intra prediction
    uint8_t ii_allowed = 0;
    for (uint8_t is_base = 0; is_base < 2; is_base++) {
        for (uint8_t transition_present = 0; transition_present < 2; transition_present++) {
            if (ii_allowed)
                break;
            ii_allowed |= svt_aom_get_inter_intra_level(enc_mode, is_base, transition_present);
        }
    }
    if (ii_allowed) {
        const uint8_t bits = ctx->hbd_md > EB_8_BIT_MD ? 2 : 1;
        // MAX block size for inter intra is 32x32
        EB_MALLOC_2D(ctx->intrapred_buf, INTERINTRA_MODES, 32 * 32 * bits * sizeof(ctx->intrapred_buf[0][0]));
    }

    // Allocate buffers for obmc prediction
    uint8_t obmc_allowed = 0;
    for (uint8_t is_ref = 0; is_ref < 2; is_ref++) {
        for (uint8_t fast_decode = 0; fast_decode < 2; fast_decode++) {
            for (EbInputResolution input_resolution = 0; input_resolution < INPUT_SIZE_COUNT; input_resolution++) {
                if (obmc_allowed)
                    break;
                obmc_allowed |= svt_aom_get_obmc_level(enc_mode, fast_decode, input_resolution);
            }
        }
    }
    if (obmc_allowed) {
        const uint8_t bits = ctx->hbd_md > EB_8_BIT_MD ? 2 : 1;
        EB_MALLOC(ctx->obmc_buff_0, sb_size * sb_size * bits * MAX_MB_PLANE * sizeof(ctx->obmc_buff_0[0]));
        EB_MALLOC(ctx->obmc_buff_1, sb_size * sb_size * bits * MAX_MB_PLANE * sizeof(ctx->obmc_buff_1[0]));
        EB_MALLOC(ctx->wsrc_buf, sb_size * sb_size * sizeof(ctx->wsrc_buf[0]));
        EB_MALLOC(ctx->mask_buf, sb_size * sb_size * sizeof(ctx->mask_buf[0]));
    }
    EB_MALLOC_ARRAY(ctx->md_local_blk_unit, block_max_count_sb);
    EB_MALLOC_ARRAY(ctx->md_blk_arr_nsq, block_max_count_sb);
    // Fast Candidate Array
    uint16_t max_can_count = svt_aom_get_max_can_count(enc_mode) + ind_uv_cands;
    EB_MALLOC_ARRAY(ctx->fast_cand_array, max_can_count);

    EB_MALLOC_ARRAY(ctx->fast_cand_ptr_array, max_can_count);
    svt_aom_assert_err(max_can_count > ind_uv_cands, "Max. candidates is too low");
    EB_MALLOC_2D(ctx->injected_mvs, (uint16_t)(max_can_count - ind_uv_cands), 2);
    EB_MALLOC_ARRAY(ctx->injected_ref_types, (max_can_count - ind_uv_cands));

    for (cand_index = 0; cand_index < max_can_count; ++cand_index) {
        ctx->fast_cand_ptr_array[cand_index]               = &ctx->fast_cand_array[cand_index];
        ctx->fast_cand_ptr_array[cand_index]->palette_info = NULL;
    }

    for (int cd = 0; cd < MAX_PAL_CAND; cd++)
        if (cfg_palette)
            EB_MALLOC_ARRAY(ctx->palette_cand_array[cd].color_idx_map, MAX_PALETTE_SQUARE);
        else
            ctx->palette_cand_array[cd].color_idx_map = NULL;

    // MD palette search
    EB_MALLOC_ARRAY(ctx->palette_size_array_0, MAX_PAL_CAND);

    // Cost Arrays
    EB_MALLOC_ARRAY(ctx->fast_cost_array, ctx->max_nics_uv);
    EB_MALLOC_ARRAY(ctx->full_cost_array, ctx->max_nics_uv);
    EB_MALLOC_ARRAY(ctx->full_cost_ssim_array, ctx->max_nics_uv);
    // Candidate Buffers
    EB_NEW(ctx->cand_bf_tx_depth_1,
           svt_aom_mode_decision_scratch_cand_bf_ctor,
           sb_size,
           ctx->hbd_md ? EB_TEN_BIT : EB_EIGHT_BIT);

    EB_ALLOC_PTR_ARRAY(ctx->cand_bf_tx_depth_1->cand, 1);
    EB_NEW(ctx->cand_bf_tx_depth_2,
           svt_aom_mode_decision_scratch_cand_bf_ctor,
           sb_size,
           ctx->hbd_md ? EB_TEN_BIT : EB_EIGHT_BIT);

    EB_ALLOC_PTR_ARRAY(ctx->cand_bf_tx_depth_2->cand, 1);
    for (int i = 0; i < 3; i++) {
        ctx->md_local_blk_unit[0].neigh_left_recon[i]       = NULL;
        ctx->md_local_blk_unit[0].neigh_top_recon[i]        = NULL;
        ctx->md_local_blk_unit[0].neigh_left_recon_16bit[i] = NULL;
        ctx->md_local_blk_unit[0].neigh_top_recon_16bit[i]  = NULL;
    }
    uint16_t sz = sizeof(uint16_t);
    if (ctx->hbd_md > EB_8_BIT_MD) {
        EB_MALLOC_ARRAY(ctx->md_local_blk_unit[0].neigh_left_recon_16bit[0], block_max_count_sb * sb_size * sz);
        EB_MALLOC_ARRAY(ctx->md_local_blk_unit[0].neigh_top_recon_16bit[0], block_max_count_sb * sb_size * sz);
        EB_MALLOC_ARRAY(ctx->md_local_blk_unit[0].neigh_left_recon_16bit[1], block_max_count_sb * sb_size * sz >> 1);
        EB_MALLOC_ARRAY(ctx->md_local_blk_unit[0].neigh_top_recon_16bit[1], block_max_count_sb * sb_size * sz >> 1);
        EB_MALLOC_ARRAY(ctx->md_local_blk_unit[0].neigh_left_recon_16bit[2], block_max_count_sb * sb_size * sz >> 1);
        EB_MALLOC_ARRAY(ctx->md_local_blk_unit[0].neigh_top_recon_16bit[2], block_max_count_sb * sb_size * sz >> 1);
    }
    if (ctx->hbd_md != EB_10_BIT_MD) {
        EB_MALLOC_ARRAY(ctx->md_local_blk_unit[0].neigh_left_recon[0], block_max_count_sb * sb_size);
        EB_MALLOC_ARRAY(ctx->md_local_blk_unit[0].neigh_top_recon[0], block_max_count_sb * sb_size);
        EB_MALLOC_ARRAY(ctx->md_local_blk_unit[0].neigh_left_recon[1], block_max_count_sb * sb_size >> 1);
        EB_MALLOC_ARRAY(ctx->md_local_blk_unit[0].neigh_top_recon[1], block_max_count_sb * sb_size >> 1);
        EB_MALLOC_ARRAY(ctx->md_local_blk_unit[0].neigh_left_recon[2], block_max_count_sb * sb_size >> 1);
        EB_MALLOC_ARRAY(ctx->md_local_blk_unit[0].neigh_top_recon[2], block_max_count_sb * sb_size >> 1);
    }
    uint32_t coded_leaf_index;
    for (coded_leaf_index = 0; coded_leaf_index < block_max_count_sb; ++coded_leaf_index) {
        size_t offset = coded_leaf_index * sb_size * sz;
        ctx->md_local_blk_unit[coded_leaf_index].neigh_left_recon_16bit[0] =
            ctx->md_local_blk_unit[0].neigh_left_recon_16bit[0] + offset;
        ctx->md_local_blk_unit[coded_leaf_index].neigh_top_recon_16bit[0] =
            ctx->md_local_blk_unit[0].neigh_top_recon_16bit[0] + offset;
        offset >>= 1;
        ctx->md_local_blk_unit[coded_leaf_index].neigh_left_recon_16bit[1] =
            ctx->md_local_blk_unit[0].neigh_left_recon_16bit[1] + offset;
        ctx->md_local_blk_unit[coded_leaf_index].neigh_top_recon_16bit[1] =
            ctx->md_local_blk_unit[0].neigh_top_recon_16bit[1] + offset;
        ctx->md_local_blk_unit[coded_leaf_index].neigh_left_recon_16bit[2] =
            ctx->md_local_blk_unit[0].neigh_left_recon_16bit[2] + offset;
        ctx->md_local_blk_unit[coded_leaf_index].neigh_top_recon_16bit[2] =
            ctx->md_local_blk_unit[0].neigh_top_recon_16bit[2] + offset;

        offset                                                       = coded_leaf_index * sb_size;
        ctx->md_local_blk_unit[coded_leaf_index].neigh_left_recon[0] = ctx->md_local_blk_unit[0].neigh_left_recon[0] +
            offset;
        ctx->md_local_blk_unit[coded_leaf_index].neigh_top_recon[0] = ctx->md_local_blk_unit[0].neigh_top_recon[0] +
            offset;
        offset >>= 1;
        ctx->md_local_blk_unit[coded_leaf_index].neigh_left_recon[1] = ctx->md_local_blk_unit[0].neigh_left_recon[1] +
            offset;
        ctx->md_local_blk_unit[coded_leaf_index].neigh_top_recon[1] = ctx->md_local_blk_unit[0].neigh_top_recon[1] +
            offset;
        ctx->md_local_blk_unit[coded_leaf_index].neigh_left_recon[2] = ctx->md_local_blk_unit[0].neigh_left_recon[2] +
            offset;
        ctx->md_local_blk_unit[coded_leaf_index].neigh_top_recon[2] = ctx->md_local_blk_unit[0].neigh_top_recon[2] +
            offset;
    }
    ctx->md_blk_arr_nsq[0].av1xd = NULL;
    EB_MALLOC_ARRAY(ctx->md_blk_arr_nsq[0].av1xd, block_max_count_sb);
    EB_MALLOC_ARRAY(ctx->avail_blk_flag, block_max_count_sb);
    EB_MALLOC_ARRAY(ctx->cost_avail, block_max_count_sb);
    EB_MALLOC_ARRAY(ctx->mdc_sb_array, 1);
    for (coded_leaf_index = 0; coded_leaf_index < block_max_count_sb; ++coded_leaf_index) {
        ctx->md_blk_arr_nsq[coded_leaf_index].av1xd      = ctx->md_blk_arr_nsq[0].av1xd + coded_leaf_index;
        ctx->md_blk_arr_nsq[coded_leaf_index].segment_id = 0;
        const BlockGeom *blk_geom                        = get_blk_geom_mds(coded_leaf_index);

        if (svt_aom_get_bypass_encdec(enc_mode, ctx->hbd_md, encoder_bit_depth)) {
            EbPictureBufferDescInitData init_data;

            init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
            init_data.max_width          = blk_geom->bwidth;
            init_data.max_height         = blk_geom->bheight;
            init_data.bit_depth          = EB_THIRTYTWO_BIT;
            init_data.color_format       = (blk_geom->bwidth > 4 && blk_geom->bheight > 4)
                      ? EB_YUV420
                      : EB_YUV444; // PW - must have at least 4x4 for chroma coeffs
            init_data.left_padding       = 0;
            init_data.right_padding      = 0;
            init_data.top_padding        = 0;
            init_data.bot_padding        = 0;
            init_data.split_mode         = FALSE;

            EB_NEW(ctx->md_blk_arr_nsq[coded_leaf_index].coeff_tmp, svt_picture_buffer_desc_ctor, (EbPtr)&init_data);

            init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
            init_data.max_width          = blk_geom->bwidth;
            init_data.max_height         = blk_geom->bheight;
            init_data.bit_depth          = ctx->hbd_md ? EB_TEN_BIT : EB_EIGHT_BIT;
            ;
            init_data.color_format  = (blk_geom->bwidth > 4 && blk_geom->bheight > 4) ? EB_YUV420 : EB_YUV444;
            init_data.left_padding  = 0;
            init_data.right_padding = 0;
            init_data.top_padding   = 0;
            init_data.bot_padding   = 0;
            init_data.split_mode    = FALSE;

            EB_NEW(ctx->md_blk_arr_nsq[coded_leaf_index].recon_tmp, svt_picture_buffer_desc_ctor, (EbPtr)&init_data);
        } else {
            ctx->md_blk_arr_nsq[coded_leaf_index].coeff_tmp = NULL;
            ctx->md_blk_arr_nsq[coded_leaf_index].recon_tmp = NULL;
        }
    }
    for (CandClass cand_class_it = CAND_CLASS_0; cand_class_it < CAND_CLASS_TOTAL; cand_class_it++)
        EB_MALLOC_ARRAY(ctx->cand_buff_indices[cand_class_it], ctx->max_nics_uv);

    EB_MALLOC_ARRAY(ctx->best_candidate_index_array, ctx->max_nics_uv);
    EB_MALLOC_ARRAY(ctx->above_txfm_context, (sb_size >> MI_SIZE_LOG2));
    EB_MALLOC_ARRAY(ctx->left_txfm_context, (sb_size >> MI_SIZE_LOG2));
    EbPictureBufferDescInitData thirty_two_width_picture_buffer_desc_init_data;
    EbPictureBufferDescInitData picture_buffer_desc_init_data;

    picture_buffer_desc_init_data.max_width          = sb_size;
    picture_buffer_desc_init_data.max_height         = sb_size;
    picture_buffer_desc_init_data.bit_depth          = ctx->hbd_md ? EB_TEN_BIT : EB_EIGHT_BIT;
    picture_buffer_desc_init_data.color_format       = EB_YUV420;
    picture_buffer_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    picture_buffer_desc_init_data.left_padding       = 0;
    picture_buffer_desc_init_data.right_padding      = 0;
    picture_buffer_desc_init_data.top_padding        = 0;
    picture_buffer_desc_init_data.bot_padding        = 0;
    picture_buffer_desc_init_data.split_mode         = FALSE;

    thirty_two_width_picture_buffer_desc_init_data.max_width          = sb_size;
    thirty_two_width_picture_buffer_desc_init_data.max_height         = sb_size;
    thirty_two_width_picture_buffer_desc_init_data.bit_depth          = EB_THIRTYTWO_BIT;
    thirty_two_width_picture_buffer_desc_init_data.color_format       = EB_YUV420;
    thirty_two_width_picture_buffer_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    thirty_two_width_picture_buffer_desc_init_data.left_padding       = 0;
    thirty_two_width_picture_buffer_desc_init_data.right_padding      = 0;
    thirty_two_width_picture_buffer_desc_init_data.top_padding        = 0;
    thirty_two_width_picture_buffer_desc_init_data.bot_padding        = 0;
    thirty_two_width_picture_buffer_desc_init_data.split_mode         = FALSE;

    for (uint32_t txt_itr = 0; txt_itr < TX_TYPES; ++txt_itr) {
        EB_NEW(ctx->recon_coeff_ptr[txt_itr],
               svt_picture_buffer_desc_ctor,
               (EbPtr)&thirty_two_width_picture_buffer_desc_init_data);
        EB_NEW(ctx->recon_ptr[txt_itr], svt_picture_buffer_desc_ctor, (EbPtr)&picture_buffer_desc_init_data);
        EB_NEW(ctx->quant_coeff_ptr[txt_itr],
               svt_picture_buffer_desc_ctor,
               (EbPtr)&thirty_two_width_picture_buffer_desc_init_data);
    }
    EB_NEW(ctx->tx_coeffs, svt_picture_buffer_desc_ctor, (EbPtr)&thirty_two_width_picture_buffer_desc_init_data);
    EB_NEW(ctx->scratch_prediction_ptr, svt_picture_buffer_desc_ctor, (EbPtr)&picture_buffer_desc_init_data);
    EbPictureBufferDescInitData double_width_picture_buffer_desc_init_data;
    double_width_picture_buffer_desc_init_data.max_width          = sb_size;
    double_width_picture_buffer_desc_init_data.max_height         = sb_size;
    double_width_picture_buffer_desc_init_data.bit_depth          = EB_SIXTEEN_BIT;
    double_width_picture_buffer_desc_init_data.color_format       = EB_YUV420;
    double_width_picture_buffer_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    double_width_picture_buffer_desc_init_data.left_padding       = 0;
    double_width_picture_buffer_desc_init_data.right_padding      = 0;
    double_width_picture_buffer_desc_init_data.top_padding        = 0;
    double_width_picture_buffer_desc_init_data.bot_padding        = 0;
    double_width_picture_buffer_desc_init_data.split_mode         = FALSE;

    // The temp_recon_ptr and temp_residual will be shared by all candidates
    // If you want to do something with residual or recon, you need to create one
    EB_NEW(ctx->temp_recon_ptr, svt_picture_buffer_desc_ctor, (EbPtr)&picture_buffer_desc_init_data);
    EB_NEW(ctx->temp_residual, svt_picture_buffer_desc_ctor, (EbPtr)&double_width_picture_buffer_desc_init_data);

    // Candidate Buffers
    EB_ALLOC_PTR_ARRAY(ctx->cand_bf_ptr_array, ctx->max_nics_uv);

    for (buffer_index = 0; buffer_index < ctx->max_nics; ++buffer_index) {
        EB_NEW(ctx->cand_bf_ptr_array[buffer_index],
               svt_aom_mode_decision_cand_bf_ctor,
               ctx->hbd_md ? EB_TEN_BIT : EB_EIGHT_BIT,
               sb_size,
               PICTURE_BUFFER_DESC_FULL_MASK,
               ctx->temp_residual,
               ctx->temp_recon_ptr,
               &(ctx->fast_cost_array[buffer_index]),
               &(ctx->full_cost_array[buffer_index]),
               &(ctx->full_cost_ssim_array[buffer_index]));
    }

    for (buffer_index = max_nics; buffer_index < ctx->max_nics_uv; ++buffer_index) {
        EB_NEW(ctx->cand_bf_ptr_array[buffer_index],
               svt_aom_mode_decision_cand_bf_ctor,
               ctx->hbd_md ? EB_TEN_BIT : EB_EIGHT_BIT,
               sb_size,
               PICTURE_BUFFER_DESC_CHROMA_MASK,
               ctx->temp_residual,
               ctx->temp_recon_ptr,
               &(ctx->fast_cost_array[buffer_index]),
               &(ctx->full_cost_array[buffer_index]),
               &(ctx->full_cost_ssim_array[buffer_index]));
    }

    return EB_ErrorNone;
}

/**************************************************
 * Reset Mode Decision Neighbor Arrays
 *************************************************/
void svt_aom_reset_mode_decision_neighbor_arrays(PictureControlSet *pcs, uint16_t tile_idx) {
    uint8_t depth;
    for (depth = 0; depth < NA_TOT_CNT; depth++) {
        svt_aom_neighbor_array_unit_reset(pcs->mdleaf_partition_na[depth][tile_idx]);
        if (pcs->hbd_md != EB_10_BIT_MD) {
            svt_aom_neighbor_array_unit_reset(pcs->md_luma_recon_na[depth][tile_idx]);
            svt_aom_neighbor_array_unit_reset(pcs->md_tx_depth_1_luma_recon_na[depth][tile_idx]);
            svt_aom_neighbor_array_unit_reset(pcs->md_tx_depth_2_luma_recon_na[depth][tile_idx]);
            svt_aom_neighbor_array_unit_reset(pcs->md_cb_recon_na[depth][tile_idx]);
            svt_aom_neighbor_array_unit_reset(pcs->md_cr_recon_na[depth][tile_idx]);
        }
        if (pcs->hbd_md > EB_8_BIT_MD) {
            svt_aom_neighbor_array_unit_reset(pcs->md_luma_recon_na_16bit[depth][tile_idx]);
            svt_aom_neighbor_array_unit_reset(pcs->md_tx_depth_1_luma_recon_na_16bit[depth][tile_idx]);
            svt_aom_neighbor_array_unit_reset(pcs->md_tx_depth_2_luma_recon_na_16bit[depth][tile_idx]);
            svt_aom_neighbor_array_unit_reset(pcs->md_cb_recon_na_16bit[depth][tile_idx]);
            svt_aom_neighbor_array_unit_reset(pcs->md_cr_recon_na_16bit[depth][tile_idx]);
        }

        svt_aom_neighbor_array_unit_reset(pcs->md_y_dcs_na[depth][tile_idx]);
        svt_aom_neighbor_array_unit_reset(pcs->md_tx_depth_1_luma_dc_sign_level_coeff_na[depth][tile_idx]);
        svt_aom_neighbor_array_unit_reset(pcs->md_cb_dc_sign_level_coeff_na[depth][tile_idx]);
        svt_aom_neighbor_array_unit_reset(pcs->md_cr_dc_sign_level_coeff_na[depth][tile_idx]);
        svt_aom_neighbor_array_unit_reset(pcs->md_txfm_context_array[depth][tile_idx]);
    }

    return;
}
// If the ref intra percentage is below the TH, applying modulation to the MD lambda
#define LAMBDA_MOD_INTRA_TH 65
#define LAMBDA_MOD_SCALING_FACTOR 138
// Set the lambda for each sb.
// When lambda tuning is on (blk_lambda_tuning), lambda of each block is set separately (full_lambda_md/fast_lambda_md)
// later in svt_aom_set_tuned_blk_lambda
// Testing showed that updating SAD lambda based on frame info was not helpful; therefore, the SAD lambda generation is not changed.
static void av1_lambda_assign_md(PictureControlSet *pcs, ModeDecisionContext *ctx) {
    ctx->full_lambda_md[0] = (uint32_t)svt_aom_compute_rd_mult(pcs, ctx->qp_index, ctx->me_q_index, 8);
    ctx->fast_lambda_md[0] = av1_lambda_mode_decision8_bit_sad[ctx->qp_index];
    ctx->full_lambda_md[1] = (uint32_t)svt_aom_compute_rd_mult(pcs, ctx->qp_index, ctx->me_q_index, 10);
    ctx->fast_lambda_md[1] = av1lambda_mode_decision10_bit_sad[ctx->qp_index];

    if (pcs->scs->stats_based_sb_lambda_modulation) {
        if (pcs->temporal_layer_index > 0) {
            if (pcs->ref_intra_percentage < LAMBDA_MOD_INTRA_TH) {
                ctx->full_lambda_md[0] = (ctx->full_lambda_md[0] * LAMBDA_MOD_SCALING_FACTOR) >> 7;
                ctx->fast_lambda_md[0] = (ctx->fast_lambda_md[0] * LAMBDA_MOD_SCALING_FACTOR) >> 7;
                ctx->full_lambda_md[1] = (ctx->full_lambda_md[1] * LAMBDA_MOD_SCALING_FACTOR) >> 7;
                ctx->fast_lambda_md[1] = (ctx->fast_lambda_md[1] * LAMBDA_MOD_SCALING_FACTOR) >> 7;
            }
        }
    }

    ctx->full_lambda_md[1] *= 16;
    ctx->fast_lambda_md[1] *= 4;

    SequenceControlSet *scs          = pcs->scs;
    uint64_t            scale_factor = scs->static_config.lambda_scale_factors[pcs->ppcs->update_type];
    ctx->full_lambda_md[0]           = (uint32_t)((ctx->full_lambda_md[0] * scale_factor) >> 7);
    ctx->full_lambda_md[1]           = (uint32_t)((ctx->full_lambda_md[1] * scale_factor) >> 7);
    ctx->fast_lambda_md[0]           = (uint32_t)((ctx->fast_lambda_md[0] * scale_factor) >> 7);
    ctx->fast_lambda_md[1]           = (uint32_t)((ctx->fast_lambda_md[1] * scale_factor) >> 7);

    ctx->full_sb_lambda_md[0] = ctx->full_lambda_md[0];
    ctx->full_sb_lambda_md[1] = ctx->full_lambda_md[1];
}

static void av1_lambda_assign(PictureControlSet *pcs, uint32_t *fast_lambda, uint32_t *full_lambda, uint8_t bit_depth,
                              uint16_t qp_index, Bool multiply_lambda) {
    if (bit_depth == 8) {
        *full_lambda = (uint32_t)svt_aom_compute_rd_mult(pcs, (uint8_t)qp_index, (uint8_t)qp_index, bit_depth);
        *fast_lambda = av1_lambda_mode_decision8_bit_sad[qp_index];
    } else if (bit_depth == 10) {
        *full_lambda = (uint32_t)svt_aom_compute_rd_mult(pcs, (uint8_t)qp_index, (uint8_t)qp_index, bit_depth);
        *fast_lambda = av1lambda_mode_decision10_bit_sad[qp_index];
        if (multiply_lambda) {
            *full_lambda *= 16;
            *fast_lambda *= 4;
        }
    } else if (bit_depth == 12) {
        *full_lambda = (uint32_t)svt_aom_compute_rd_mult(pcs, (uint8_t)qp_index, (uint8_t)qp_index, bit_depth);
        *fast_lambda = av1lambda_mode_decision12_bit_sad[qp_index];
    } else {
        assert(bit_depth >= 8);
        assert(bit_depth <= 12);
    }

    // NM: To be done: tune lambda based on the picture type and layer.
    SequenceControlSet *scs          = pcs->scs;
    uint64_t            scale_factor = scs->static_config.lambda_scale_factors[pcs->ppcs->update_type];
    *full_lambda                     = (uint32_t)((*full_lambda * scale_factor) >> 7);
    *fast_lambda                     = (uint32_t)((*fast_lambda * scale_factor) >> 7);
}

const EbAv1LambdaAssignFunc svt_aom_av1_lambda_assignment_function_table[4] = {
    av1_lambda_assign,
    av1_lambda_assign,
    av1_lambda_assign,
    av1_lambda_assign,
};

void svt_aom_reset_mode_decision(SequenceControlSet *scs, ModeDecisionContext *ctx, PictureControlSet *pcs,
                                 uint16_t tile_group_idx, uint32_t segment_index) {
    ctx->hbd_md = pcs->hbd_md;
    // Reset MD rate Estimation table to initial values by copying from md_rate_est_ctx
    ctx->md_rate_est_ctx = pcs->md_rate_est_ctx;
    // Reset CABAC Contexts

    // Reset Neighbor Arrays at start of new Segment / Picture
    if (segment_index == 0) {
        for (uint16_t r = pcs->ppcs->tile_group_info[tile_group_idx].tile_group_tile_start_y;
             r < pcs->ppcs->tile_group_info[tile_group_idx].tile_group_tile_end_y;
             r++) {
            for (uint16_t c = pcs->ppcs->tile_group_info[tile_group_idx].tile_group_tile_start_x;
                 c < pcs->ppcs->tile_group_info[tile_group_idx].tile_group_tile_end_x;
                 c++) {
                uint16_t tile_idx = c + r * pcs->ppcs->av1_cm->tiles_info.tile_cols;
                svt_aom_reset_mode_decision_neighbor_arrays(pcs, tile_idx);
            }
        }
        (void)scs;
    }
    //each segment enherits the bypass encdec from the picture level
    ctx->bypass_encdec = pcs->pic_bypass_encdec;
    ctx->skip_pd0      = pcs->pic_skip_pd0;
    set_block_based_depth_refinement_controls(ctx, pcs->pic_block_based_depth_refinement_level);
    if (!pcs->rtc_tune || pcs->temporal_layer_index != 0)
        ctx->rtc_use_N4_dct_dct_shortcut = 1;
    else
        ctx->rtc_use_N4_dct_dct_shortcut = 0;
    return;
}

/******************************************************
 * Mode Decision Configure SB
 ******************************************************/
void svt_aom_mode_decision_configure_sb(ModeDecisionContext *ctx, PictureControlSet *pcs, uint8_t sb_qp,
                                        uint8_t me_sb_qp) {
    /* Note(CHKN) : when Qp modulation varies QP on a sub-SB(CU) basis,  Lamda has to change based on Cu->QP , and then this code has to move inside the CU loop in MD */

    // Lambda Assignement
    ctx->qp_index = pcs->ppcs->frm_hdr.delta_q_params.delta_q_present
        ? sb_qp
        : (uint8_t)pcs->ppcs->frm_hdr.quantization_params.base_q_idx;

    ctx->me_q_index = me_sb_qp;

    av1_lambda_assign_md(pcs, ctx);

    ctx->hbd_pack_done = 0;

    return;
}
