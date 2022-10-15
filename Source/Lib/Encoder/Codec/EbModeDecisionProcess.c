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

uint8_t get_bypass_encdec(EncMode enc_mode, uint8_t hbd_mode_decision, uint8_t encoder_bit_depth);
void    set_block_based_depth_refinement_controls(ModeDecisionContext *mdctxt,
                                                  uint8_t block_based_depth_refinement_level);
int     svt_av1_allow_palette(int allow_palette, BlockSize sb_type);
static void mode_decision_context_dctor(EbPtr p) {
    ModeDecisionContext *obj = (ModeDecisionContext *)p;

    uint32_t block_max_count_sb = obj->init_max_block_cnt;
    for (int cd = 0; cd < MAX_PAL_CAND; cd++)
        if (obj->palette_cand_array[cd].color_idx_map)
            EB_FREE_ARRAY(obj->palette_cand_array[cd].color_idx_map);

    // MD palette search
    if (obj->palette_size_array_0)
        EB_FREE_ARRAY(obj->palette_size_array_0);
    if (obj->palette_size_array_1)
        EB_FREE_ARRAY(obj->palette_size_array_1);
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
    EB_DELETE_PTR_ARRAY(obj->candidate_buffer_ptr_array, obj->max_nics_uv);
    EB_FREE_ARRAY(obj->candidate_buffer_tx_depth_1->candidate_ptr);
    EB_DELETE(obj->candidate_buffer_tx_depth_1);
    EB_FREE_ARRAY(obj->candidate_buffer_tx_depth_2->candidate_ptr);
    EB_DELETE(obj->candidate_buffer_tx_depth_2);
    EB_FREE_ALIGNED_ARRAY(obj->cfl_temp_luma_recon16bit);
    EB_FREE_ALIGNED_ARRAY(obj->cfl_temp_luma_recon);
    EB_FREE_ARRAY(obj->fast_candidate_array);
    EB_FREE_ARRAY(obj->fast_candidate_ptr_array);
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
    EB_FREE_ARRAY(obj->tested_blk_flag);
    EB_FREE_ARRAY(obj->do_not_process_blk);
    EB_FREE_ARRAY(obj->md_local_blk_unit);
    EB_FREE_ARRAY(obj->md_blk_arr_nsq);
    if (obj->rate_est_table)
        EB_FREE_ARRAY(obj->rate_est_table);
    EB_FREE_ARRAY(obj->mdc_sb_array);
    for (uint32_t txt_itr = 0; txt_itr < TX_TYPES; ++txt_itr) {
        EB_DELETE(obj->recon_coeff_ptr[txt_itr]);
        EB_DELETE(obj->recon_ptr[txt_itr]);
        EB_DELETE(obj->quant_coeff_ptr[txt_itr]);
    }
    EB_DELETE(obj->tx_coeffs);
    EB_DELETE(obj->scratch_prediction_ptr);
    EB_DELETE(obj->temp_residual_ptr);
    EB_DELETE(obj->temp_recon_ptr);
}
uint8_t get_nic_level(EncMode enc_mode, uint8_t is_base, uint8_t hierarchical_levels);
uint8_t set_nic_controls(ModeDecisionContext *ctx, uint8_t nic_level);
void    set_nics(NicScalingCtrls *scaling_ctrls, uint32_t mds1_count[CAND_CLASS_TOTAL],
                 uint32_t mds2_count[CAND_CLASS_TOTAL], uint32_t mds3_count[CAND_CLASS_TOTAL],
                 uint8_t pic_type);

uint8_t get_update_cdf_level(EncMode enc_mode, SliceType is_islice, uint8_t is_base);
uint8_t svt_aom_get_chroma_level(EncMode enc_mode);
uint8_t svt_aom_set_chroma_controls(ModeDecisionContext *ctx, uint8_t uv_level);

/*
* return the max canidate count for MDS0
  Used by candidate injection and memory allocation
*/
uint16_t get_max_can_count(EncMode enc_mode) {
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
 * Mode Decision Context Constructor
 ******************************************************/
EbErrorType mode_decision_context_ctor(ModeDecisionContext *context_ptr, EbColorFormat color_format,
                                       uint8_t sb_size, uint8_t enc_mode, uint16_t max_block_cnt,
                                       uint32_t encoder_bit_depth,
                                       EbFifo  *mode_decision_configuration_input_fifo_ptr,
                                       EbFifo  *mode_decision_output_fifo_ptr,
                                       uint8_t enable_hbd_mode_decision, uint8_t cfg_palette,
                                       uint32_t hierarchical_levels) {
    uint32_t buffer_index;
    uint32_t cand_index;

    context_ptr->init_max_block_cnt = max_block_cnt;
    uint32_t block_max_count_sb     = max_block_cnt;

    context_ptr->sb_size = sb_size;
    (void)color_format;

    context_ptr->dctor             = mode_decision_context_dctor;
    context_ptr->hbd_mode_decision = enable_hbd_mode_decision;

    // Input/Output System Resource Manager FIFOs
    context_ptr->mode_decision_configuration_input_fifo_ptr =
        mode_decision_configuration_input_fifo_ptr;
    context_ptr->mode_decision_output_fifo_ptr = mode_decision_output_fifo_ptr;

    // Maximum number of candidates MD can support
    // determine MAX_NICS for a given preset
    // get the min scaling level (the smallest scaling level is the most conservative)
    uint8_t min_nic_scaling_level = NICS_SCALING_LEVELS - 1;
    for (uint8_t is_base = 0; is_base < 2; is_base++) {
        uint8_t nic_level         = get_nic_level(enc_mode, is_base, hierarchical_levels);
        uint8_t nic_scaling_level = set_nic_controls(NULL, nic_level);
        min_nic_scaling_level     = MIN(min_nic_scaling_level, nic_scaling_level);
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
            set_nics(&scaling_ctrls, mds1_count, mds2_count, mds3_count, pic_type);

            uint32_t nics = 0;
            for (CandClass cidx = CAND_CLASS_0; cidx < CAND_CLASS_TOTAL; cidx++) {
                nics += mds1_count[cidx];
            }
            max_nics = MAX(max_nics, nics);
        }
    }

    // If independent chroma search is used, need to allocate additional 84 candidate buffers
    const uint8_t ind_uv_cands = svt_aom_set_chroma_controls(
                                     NULL, svt_aom_get_chroma_level(enc_mode)) == CHROMA_MODE_0
        ? 84
        : 0;
    max_nics += CAND_CLASS_TOTAL; //need one extra temp buffer for each fast loop call
    context_ptr->max_nics    = max_nics;
    context_ptr->max_nics_uv = max_nics + ind_uv_cands;
    // Cfl scratch memory
    if (context_ptr->hbd_mode_decision > EB_8_BIT_MD)
        EB_MALLOC_ALIGNED(context_ptr->cfl_temp_luma_recon16bit,
                          sizeof(uint16_t) * sb_size * sb_size);
    if (context_ptr->hbd_mode_decision != EB_10_BIT_MD)
        EB_MALLOC_ALIGNED(context_ptr->cfl_temp_luma_recon, sizeof(uint8_t) * sb_size * sb_size);
    uint8_t use_update_cdf = 0;
    for (uint8_t is_islice = 0; is_islice < 2; is_islice++) {
        for (uint8_t is_base = 0; is_base < 2; is_base++) {
            if (use_update_cdf)
                break;
            use_update_cdf |= get_update_cdf_level(enc_mode, is_islice, is_base);
        }
    }
    if (use_update_cdf)
        EB_CALLOC_ARRAY(context_ptr->rate_est_table, 1);
    else
        context_ptr->rate_est_table = NULL;
    EB_MALLOC_ARRAY(context_ptr->md_local_blk_unit, block_max_count_sb);
    EB_MALLOC_ARRAY(context_ptr->md_blk_arr_nsq, block_max_count_sb);
    // Fast Candidate Array
    uint16_t max_can_count = get_max_can_count(enc_mode) + ind_uv_cands;
    EB_MALLOC_ARRAY(context_ptr->fast_candidate_array, max_can_count);

    EB_MALLOC_ARRAY(context_ptr->fast_candidate_ptr_array, max_can_count);
    assert_err(max_can_count > ind_uv_cands, "Max. candidates is too low");
    EB_MALLOC_2D(context_ptr->injected_mvs, (uint16_t)(max_can_count - ind_uv_cands), 2);
    EB_MALLOC_ARRAY(context_ptr->injected_ref_types, (max_can_count - ind_uv_cands));

    for (cand_index = 0; cand_index < max_can_count; ++cand_index) {
        context_ptr->fast_candidate_ptr_array[cand_index] =
            &context_ptr->fast_candidate_array[cand_index];
        context_ptr->fast_candidate_ptr_array[cand_index]->palette_info = NULL;
    }

    for (int cd = 0; cd < MAX_PAL_CAND; cd++)
        if (cfg_palette)
            EB_MALLOC_ARRAY(context_ptr->palette_cand_array[cd].color_idx_map, MAX_PALETTE_SQUARE);
        else
            context_ptr->palette_cand_array[cd].color_idx_map = NULL;

    // MD palette search
    EB_MALLOC_ARRAY(context_ptr->palette_size_array_0, MAX_PAL_CAND);
    EB_MALLOC_ARRAY(context_ptr->palette_size_array_1, MAX_PAL_CAND);

    // Cost Arrays
    EB_MALLOC_ARRAY(context_ptr->fast_cost_array, context_ptr->max_nics_uv);
    EB_MALLOC_ARRAY(context_ptr->full_cost_array, context_ptr->max_nics_uv);
    // Candidate Buffers
    EB_NEW(context_ptr->candidate_buffer_tx_depth_1,
           mode_decision_scratch_candidate_buffer_ctor,
           sb_size,
           context_ptr->hbd_mode_decision ? EB_TEN_BIT : EB_EIGHT_BIT);

    EB_ALLOC_PTR_ARRAY(context_ptr->candidate_buffer_tx_depth_1->candidate_ptr, 1);
    EB_NEW(context_ptr->candidate_buffer_tx_depth_2,
           mode_decision_scratch_candidate_buffer_ctor,
           sb_size,
           context_ptr->hbd_mode_decision ? EB_TEN_BIT : EB_EIGHT_BIT);

    EB_ALLOC_PTR_ARRAY(context_ptr->candidate_buffer_tx_depth_2->candidate_ptr, 1);
    for (int i = 0; i < 3; i++) {
        context_ptr->md_local_blk_unit[0].neigh_left_recon[i]       = NULL;
        context_ptr->md_local_blk_unit[0].neigh_top_recon[i]        = NULL;
        context_ptr->md_local_blk_unit[0].neigh_left_recon_16bit[i] = NULL;
        context_ptr->md_local_blk_unit[0].neigh_top_recon_16bit[i]  = NULL;
    }
    uint16_t sz = sizeof(uint16_t);
    if (context_ptr->hbd_mode_decision > EB_8_BIT_MD) {
        EB_MALLOC_ARRAY(context_ptr->md_local_blk_unit[0].neigh_left_recon_16bit[0],
                        block_max_count_sb * sb_size * sz);
        EB_MALLOC_ARRAY(context_ptr->md_local_blk_unit[0].neigh_top_recon_16bit[0],
                        block_max_count_sb * sb_size * sz);
        EB_MALLOC_ARRAY(context_ptr->md_local_blk_unit[0].neigh_left_recon_16bit[1],
                        block_max_count_sb * sb_size * sz >> 1);
        EB_MALLOC_ARRAY(context_ptr->md_local_blk_unit[0].neigh_top_recon_16bit[1],
                        block_max_count_sb * sb_size * sz >> 1);
        EB_MALLOC_ARRAY(context_ptr->md_local_blk_unit[0].neigh_left_recon_16bit[2],
                        block_max_count_sb * sb_size * sz >> 1);
        EB_MALLOC_ARRAY(context_ptr->md_local_blk_unit[0].neigh_top_recon_16bit[2],
                        block_max_count_sb * sb_size * sz >> 1);
    }
    if (context_ptr->hbd_mode_decision != EB_10_BIT_MD) {
        EB_MALLOC_ARRAY(context_ptr->md_local_blk_unit[0].neigh_left_recon[0],
                        block_max_count_sb * sb_size);
        EB_MALLOC_ARRAY(context_ptr->md_local_blk_unit[0].neigh_top_recon[0],
                        block_max_count_sb * sb_size);
        EB_MALLOC_ARRAY(context_ptr->md_local_blk_unit[0].neigh_left_recon[1],
                        block_max_count_sb * sb_size >> 1);
        EB_MALLOC_ARRAY(context_ptr->md_local_blk_unit[0].neigh_top_recon[1],
                        block_max_count_sb * sb_size >> 1);
        EB_MALLOC_ARRAY(context_ptr->md_local_blk_unit[0].neigh_left_recon[2],
                        block_max_count_sb * sb_size >> 1);
        EB_MALLOC_ARRAY(context_ptr->md_local_blk_unit[0].neigh_top_recon[2],
                        block_max_count_sb * sb_size >> 1);
    }
    uint32_t coded_leaf_index;
    for (coded_leaf_index = 0; coded_leaf_index < block_max_count_sb; ++coded_leaf_index) {
        size_t offset = coded_leaf_index * sb_size * sz;
        context_ptr->md_local_blk_unit[coded_leaf_index].neigh_left_recon_16bit[0] =
            context_ptr->md_local_blk_unit[0].neigh_left_recon_16bit[0] + offset;
        context_ptr->md_local_blk_unit[coded_leaf_index].neigh_top_recon_16bit[0] =
            context_ptr->md_local_blk_unit[0].neigh_top_recon_16bit[0] + offset;
        offset >>= 1;
        context_ptr->md_local_blk_unit[coded_leaf_index].neigh_left_recon_16bit[1] =
            context_ptr->md_local_blk_unit[0].neigh_left_recon_16bit[1] + offset;
        context_ptr->md_local_blk_unit[coded_leaf_index].neigh_top_recon_16bit[1] =
            context_ptr->md_local_blk_unit[0].neigh_top_recon_16bit[1] + offset;
        context_ptr->md_local_blk_unit[coded_leaf_index].neigh_left_recon_16bit[2] =
            context_ptr->md_local_blk_unit[0].neigh_left_recon_16bit[2] + offset;
        context_ptr->md_local_blk_unit[coded_leaf_index].neigh_top_recon_16bit[2] =
            context_ptr->md_local_blk_unit[0].neigh_top_recon_16bit[2] + offset;

        offset = coded_leaf_index * sb_size;
        context_ptr->md_local_blk_unit[coded_leaf_index].neigh_left_recon[0] =
            context_ptr->md_local_blk_unit[0].neigh_left_recon[0] + offset;
        context_ptr->md_local_blk_unit[coded_leaf_index].neigh_top_recon[0] =
            context_ptr->md_local_blk_unit[0].neigh_top_recon[0] + offset;
        offset >>= 1;
        context_ptr->md_local_blk_unit[coded_leaf_index].neigh_left_recon[1] =
            context_ptr->md_local_blk_unit[0].neigh_left_recon[1] + offset;
        context_ptr->md_local_blk_unit[coded_leaf_index].neigh_top_recon[1] =
            context_ptr->md_local_blk_unit[0].neigh_top_recon[1] + offset;
        context_ptr->md_local_blk_unit[coded_leaf_index].neigh_left_recon[2] =
            context_ptr->md_local_blk_unit[0].neigh_left_recon[2] + offset;
        context_ptr->md_local_blk_unit[coded_leaf_index].neigh_top_recon[2] =
            context_ptr->md_local_blk_unit[0].neigh_top_recon[2] + offset;
    }
    context_ptr->md_blk_arr_nsq[0].av1xd = NULL;
    EB_MALLOC_ARRAY(context_ptr->md_blk_arr_nsq[0].av1xd, block_max_count_sb);
    EB_MALLOC_ARRAY(context_ptr->avail_blk_flag, block_max_count_sb);
    EB_MALLOC_ARRAY(context_ptr->tested_blk_flag, block_max_count_sb);
    EB_MALLOC_ARRAY(context_ptr->do_not_process_blk, block_max_count_sb);
    EB_MALLOC_ARRAY(context_ptr->mdc_sb_array, 1);
    for (coded_leaf_index = 0; coded_leaf_index < block_max_count_sb; ++coded_leaf_index) {
        context_ptr->md_blk_arr_nsq[coded_leaf_index].av1xd = context_ptr->md_blk_arr_nsq[0].av1xd +
            coded_leaf_index;
        context_ptr->md_blk_arr_nsq[coded_leaf_index].segment_id = 0;
        const BlockGeom *blk_geom = get_blk_geom_mds(coded_leaf_index);

        if (get_bypass_encdec(enc_mode, context_ptr->hbd_mode_decision, encoder_bit_depth)) {
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

            EB_NEW(context_ptr->md_blk_arr_nsq[coded_leaf_index].coeff_tmp,
                   svt_picture_buffer_desc_ctor,
                   (EbPtr)&init_data);

            init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
            init_data.max_width          = blk_geom->bwidth;
            init_data.max_height         = blk_geom->bheight;
            init_data.bit_depth = context_ptr->hbd_mode_decision ? EB_TEN_BIT : EB_EIGHT_BIT;
            ;
            init_data.color_format  = (blk_geom->bwidth > 4 && blk_geom->bheight > 4) ? EB_YUV420
                                                                                      : EB_YUV444;
            init_data.left_padding  = 0;
            init_data.right_padding = 0;
            init_data.top_padding   = 0;
            init_data.bot_padding   = 0;
            init_data.split_mode    = FALSE;

            EB_NEW(context_ptr->md_blk_arr_nsq[coded_leaf_index].recon_tmp,
                   svt_picture_buffer_desc_ctor,
                   (EbPtr)&init_data);
        } else {
            context_ptr->md_blk_arr_nsq[coded_leaf_index].coeff_tmp = NULL;
            context_ptr->md_blk_arr_nsq[coded_leaf_index].recon_tmp = NULL;
        }
    }
    for (CandClass cand_class_it = CAND_CLASS_0; cand_class_it < CAND_CLASS_TOTAL; cand_class_it++)
        EB_MALLOC_ARRAY(context_ptr->cand_buff_indices[cand_class_it], context_ptr->max_nics_uv);

    EB_MALLOC_ARRAY(context_ptr->best_candidate_index_array, context_ptr->max_nics_uv);
    EB_MALLOC_ARRAY(context_ptr->above_txfm_context, (sb_size >> MI_SIZE_LOG2));
    EB_MALLOC_ARRAY(context_ptr->left_txfm_context, (sb_size >> MI_SIZE_LOG2));
    EbPictureBufferDescInitData thirty_two_width_picture_buffer_desc_init_data;
    EbPictureBufferDescInitData picture_buffer_desc_init_data;

    picture_buffer_desc_init_data.max_width          = sb_size;
    picture_buffer_desc_init_data.max_height         = sb_size;
    picture_buffer_desc_init_data.bit_depth          = context_ptr->hbd_mode_decision ? EB_TEN_BIT
                                                                                      : EB_EIGHT_BIT;
    picture_buffer_desc_init_data.color_format       = EB_YUV420;
    picture_buffer_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    picture_buffer_desc_init_data.left_padding       = 0;
    picture_buffer_desc_init_data.right_padding      = 0;
    picture_buffer_desc_init_data.top_padding        = 0;
    picture_buffer_desc_init_data.bot_padding        = 0;
    picture_buffer_desc_init_data.split_mode         = FALSE;

    thirty_two_width_picture_buffer_desc_init_data.max_width    = sb_size;
    thirty_two_width_picture_buffer_desc_init_data.max_height   = sb_size;
    thirty_two_width_picture_buffer_desc_init_data.bit_depth    = EB_THIRTYTWO_BIT;
    thirty_two_width_picture_buffer_desc_init_data.color_format = EB_YUV420;
    thirty_two_width_picture_buffer_desc_init_data.buffer_enable_mask =
        PICTURE_BUFFER_DESC_FULL_MASK;
    thirty_two_width_picture_buffer_desc_init_data.left_padding  = 0;
    thirty_two_width_picture_buffer_desc_init_data.right_padding = 0;
    thirty_two_width_picture_buffer_desc_init_data.top_padding   = 0;
    thirty_two_width_picture_buffer_desc_init_data.bot_padding   = 0;
    thirty_two_width_picture_buffer_desc_init_data.split_mode    = FALSE;

    for (uint32_t txt_itr = 0; txt_itr < TX_TYPES; ++txt_itr) {
        EB_NEW(context_ptr->recon_coeff_ptr[txt_itr],
               svt_picture_buffer_desc_ctor,
               (EbPtr)&thirty_two_width_picture_buffer_desc_init_data);
        EB_NEW(context_ptr->recon_ptr[txt_itr],
               svt_picture_buffer_desc_ctor,
               (EbPtr)&picture_buffer_desc_init_data);
        EB_NEW(context_ptr->quant_coeff_ptr[txt_itr],
               svt_picture_buffer_desc_ctor,
               (EbPtr)&thirty_two_width_picture_buffer_desc_init_data);
    }
    EB_NEW(context_ptr->tx_coeffs,
           svt_picture_buffer_desc_ctor,
           (EbPtr)&thirty_two_width_picture_buffer_desc_init_data);
    EB_NEW(context_ptr->scratch_prediction_ptr,
           svt_picture_buffer_desc_ctor,
           (EbPtr)&picture_buffer_desc_init_data);
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

    // The temp_recon_ptr and temp_residual_ptr will be shared by all candidates
    // If you want to do something with residual or recon, you need to create one
    EB_NEW(context_ptr->temp_recon_ptr,
           svt_picture_buffer_desc_ctor,
           (EbPtr)&picture_buffer_desc_init_data);
    EB_NEW(context_ptr->temp_residual_ptr,
           svt_picture_buffer_desc_ctor,
           (EbPtr)&double_width_picture_buffer_desc_init_data);

    // Candidate Buffers
    EB_ALLOC_PTR_ARRAY(context_ptr->candidate_buffer_ptr_array, context_ptr->max_nics_uv);

    for (buffer_index = 0; buffer_index < context_ptr->max_nics; ++buffer_index) {
        EB_NEW(context_ptr->candidate_buffer_ptr_array[buffer_index],
               mode_decision_candidate_buffer_ctor,
               context_ptr->hbd_mode_decision ? EB_TEN_BIT : EB_EIGHT_BIT,
               sb_size,
               PICTURE_BUFFER_DESC_FULL_MASK,
               context_ptr->temp_residual_ptr,
               context_ptr->temp_recon_ptr,
               &(context_ptr->fast_cost_array[buffer_index]),
               &(context_ptr->full_cost_array[buffer_index]));
    }

    for (buffer_index = max_nics; buffer_index < context_ptr->max_nics_uv; ++buffer_index) {
        EB_NEW(context_ptr->candidate_buffer_ptr_array[buffer_index],
               mode_decision_candidate_buffer_ctor,
               context_ptr->hbd_mode_decision ? EB_TEN_BIT : EB_EIGHT_BIT,
               sb_size,
               PICTURE_BUFFER_DESC_CHROMA_MASK,
               context_ptr->temp_residual_ptr,
               context_ptr->temp_recon_ptr,
               &(context_ptr->fast_cost_array[buffer_index]),
               &(context_ptr->full_cost_array[buffer_index]));
    }
    return EB_ErrorNone;
}

/**************************************************
 * Reset Mode Decision Neighbor Arrays
 *************************************************/
void reset_mode_decision_neighbor_arrays(PictureControlSet *pcs_ptr, uint16_t tile_idx) {
    uint8_t depth;
    for (depth = 0; depth < NEIGHBOR_ARRAY_TOTAL_COUNT; depth++) {
        neighbor_array_unit_reset(pcs_ptr->md_intra_luma_mode_neighbor_array[depth][tile_idx]);
        neighbor_array_unit_reset(pcs_ptr->md_skip_flag_neighbor_array[depth][tile_idx]);
        neighbor_array_unit_reset(pcs_ptr->md_mode_type_neighbor_array[depth][tile_idx]);
        neighbor_array_unit_reset(pcs_ptr->mdleaf_partition_neighbor_array[depth][tile_idx]);
        if (pcs_ptr->hbd_mode_decision != EB_10_BIT_MD) {
            neighbor_array_unit_reset(pcs_ptr->md_luma_recon_neighbor_array[depth][tile_idx]);
            neighbor_array_unit_reset(
                pcs_ptr->md_tx_depth_1_luma_recon_neighbor_array[depth][tile_idx]);
            neighbor_array_unit_reset(
                pcs_ptr->md_tx_depth_2_luma_recon_neighbor_array[depth][tile_idx]);
            neighbor_array_unit_reset(pcs_ptr->md_cb_recon_neighbor_array[depth][tile_idx]);
            neighbor_array_unit_reset(pcs_ptr->md_cr_recon_neighbor_array[depth][tile_idx]);
        }
        if (pcs_ptr->hbd_mode_decision > EB_8_BIT_MD) {
            neighbor_array_unit_reset(pcs_ptr->md_luma_recon_neighbor_array16bit[depth][tile_idx]);
            neighbor_array_unit_reset(
                pcs_ptr->md_tx_depth_1_luma_recon_neighbor_array16bit[depth][tile_idx]);
            neighbor_array_unit_reset(
                pcs_ptr->md_tx_depth_2_luma_recon_neighbor_array16bit[depth][tile_idx]);
            neighbor_array_unit_reset(pcs_ptr->md_cb_recon_neighbor_array16bit[depth][tile_idx]);
            neighbor_array_unit_reset(pcs_ptr->md_cr_recon_neighbor_array16bit[depth][tile_idx]);
        }

        neighbor_array_unit_reset(
            pcs_ptr->md_luma_dc_sign_level_coeff_neighbor_array[depth][tile_idx]);
        neighbor_array_unit_reset(
            pcs_ptr->md_tx_depth_1_luma_dc_sign_level_coeff_neighbor_array[depth][tile_idx]);
        neighbor_array_unit_reset(
            pcs_ptr->md_cb_dc_sign_level_coeff_neighbor_array[depth][tile_idx]);
        neighbor_array_unit_reset(
            pcs_ptr->md_cr_dc_sign_level_coeff_neighbor_array[depth][tile_idx]);
        neighbor_array_unit_reset(pcs_ptr->md_txfm_context_array[depth][tile_idx]);
        neighbor_array_unit_reset(pcs_ptr->md_skip_coeff_neighbor_array[depth][tile_idx]);
        neighbor_array_unit_reset(pcs_ptr->md_ref_frame_type_neighbor_array[depth][tile_idx]);

        neighbor_array_unit_reset32(pcs_ptr->md_interpolation_type_neighbor_array[depth][tile_idx]);
    }

    return;
}

// Set the lambda for each sb.
// When lambda tuning is on (blk_lambda_tuning), lambda of each block is set separately (full_lambda_md/fast_lambda_md)
// later in set_tuned_blk_lambda
// Testing showed that updating SAD lambda based on frame info was not helpful; therefore, the SAD lambda generation is not changed.
static void av1_lambda_assign_md(PictureControlSet *pcs, ModeDecisionContext *ctx) {
    ctx->full_lambda_md[0] = (uint32_t)svt_aom_compute_rd_mult(
        pcs, ctx->qp_index, ctx->me_q_index, 8);
    ctx->fast_lambda_md[0] = av1_lambda_mode_decision8_bit_sad[ctx->qp_index];
    ctx->full_lambda_md[1] = (uint32_t)svt_aom_compute_rd_mult(
        pcs, ctx->qp_index, ctx->me_q_index, 10);
    ctx->fast_lambda_md[1] = av1lambda_mode_decision10_bit_sad[ctx->qp_index];

    ctx->full_lambda_md[1] *= 16;
    ctx->fast_lambda_md[1] *= 4;
    ctx->full_sb_lambda_md[0] = ctx->full_lambda_md[0];
    ctx->full_sb_lambda_md[1] = ctx->full_lambda_md[1];
}

void av1_lambda_assign(PictureControlSet *pcs_ptr, uint32_t *fast_lambda, uint32_t *full_lambda,
                       uint8_t bit_depth, uint16_t qp_index, Bool multiply_lambda) {
    if (bit_depth == 8) {
        *full_lambda = (uint32_t)svt_aom_compute_rd_mult(
            pcs_ptr, (uint8_t)qp_index, (uint8_t)qp_index, bit_depth);
        *fast_lambda = av1_lambda_mode_decision8_bit_sad[qp_index];
    } else if (bit_depth == 10) {
        *full_lambda = (uint32_t)svt_aom_compute_rd_mult(
            pcs_ptr, (uint8_t)qp_index, (uint8_t)qp_index, bit_depth);
        *fast_lambda = av1lambda_mode_decision10_bit_sad[qp_index];
        if (multiply_lambda) {
            *full_lambda *= 16;
            *fast_lambda *= 4;
        }
    } else if (bit_depth == 12) {
        *full_lambda = (uint32_t)svt_aom_compute_rd_mult(
            pcs_ptr, (uint8_t)qp_index, (uint8_t)qp_index, bit_depth);
        *fast_lambda = av1lambda_mode_decision12_bit_sad[qp_index];
    } else {
        assert(bit_depth >= 8);
        assert(bit_depth <= 12);
    }

    // NM: To be done: tune lambda based on the picture type and layer.
}
const EbAv1LambdaAssignFunc av1_lambda_assignment_function_table[4] = {
    av1_lambda_assign,
    av1_lambda_assign,
    av1_lambda_assign,
    av1_lambda_assign,
};

void reset_mode_decision(SequenceControlSet *scs_ptr, ModeDecisionContext *context_ptr,
                         PictureControlSet *pcs_ptr, uint16_t tile_group_idx,
                         uint32_t segment_index) {
    context_ptr->hbd_mode_decision = pcs_ptr->hbd_mode_decision;
    // Reset MD rate Estimation table to initial values by copying from md_rate_estimation_array
    context_ptr->md_rate_estimation_ptr = pcs_ptr->md_rate_estimation_array;
    // Reset CABAC Contexts

    // Reset Neighbor Arrays at start of new Segment / Picture
    if (segment_index == 0) {
        for (uint16_t r =
                 pcs_ptr->parent_pcs_ptr->tile_group_info[tile_group_idx].tile_group_tile_start_y;
             r < pcs_ptr->parent_pcs_ptr->tile_group_info[tile_group_idx].tile_group_tile_end_y;
             r++) {
            for (uint16_t c = pcs_ptr->parent_pcs_ptr->tile_group_info[tile_group_idx]
                                  .tile_group_tile_start_x;
                 c < pcs_ptr->parent_pcs_ptr->tile_group_info[tile_group_idx].tile_group_tile_end_x;
                 c++) {
                uint16_t tile_idx = c + r * pcs_ptr->parent_pcs_ptr->av1_cm->tiles_info.tile_cols;
                reset_mode_decision_neighbor_arrays(pcs_ptr, tile_idx);
            }
        }
        (void)scs_ptr;
    }
    //each segment enherits the bypass encdec from the picture level
    context_ptr->bypass_encdec = pcs_ptr->pic_bypass_encdec;
    context_ptr->skip_pd0      = pcs_ptr->pic_skip_pd0;
    set_block_based_depth_refinement_controls(context_ptr,
                                              pcs_ptr->pic_block_based_depth_refinement_level);
    return;
}

/******************************************************
 * Mode Decision Configure SB
 ******************************************************/
void mode_decision_configure_sb(ModeDecisionContext *context_ptr, PictureControlSet *pcs_ptr,
                                uint8_t sb_qp, uint8_t me_sb_qp) {
    /* Note(CHKN) : when Qp modulation varies QP on a sub-SB(CU) basis,  Lamda has to change based on Cu->QP , and then this code has to move inside the CU loop in MD */

    // Lambda Assignement
    context_ptr->qp_index = pcs_ptr->parent_pcs_ptr->frm_hdr.delta_q_params.delta_q_present
        ? sb_qp
        : (uint8_t)pcs_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx;

    context_ptr->me_q_index = me_sb_qp;

    av1_lambda_assign_md(pcs_ptr, context_ptr);

    context_ptr->hbd_pack_done = 0;

    return;
}
