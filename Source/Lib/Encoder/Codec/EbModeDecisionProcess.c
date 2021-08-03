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

int         svt_av1_allow_palette(int allow_palette, BlockSize sb_type);
static void mode_decision_context_dctor(EbPtr p) {
    ModeDecisionContext *obj                = (ModeDecisionContext *)p;

#if CLN_GEOM
    uint32_t block_max_count_sb = obj->init_max_block_cnt;
#else
    uint32_t             block_max_count_sb = (obj->sb_size == MAX_SB_SIZE) ? BLOCK_MAX_COUNT_SB_128
                                                                            : BLOCK_MAX_COUNT_SB_64;
#endif
    for (int cd = 0; cd < MAX_PAL_CAND; cd++)
        if (obj->palette_cand_array[cd].color_idx_map)
            EB_FREE_ARRAY(obj->palette_cand_array[cd].color_idx_map);
    for (uint32_t cand_index = 0; cand_index < MODE_DECISION_CANDIDATE_MAX_COUNT; ++cand_index) {
        for (uint32_t coded_leaf_index = 0; coded_leaf_index < block_max_count_sb;
             ++coded_leaf_index)
            if (obj->md_blk_arr_nsq[coded_leaf_index].palette_info.color_idx_map)
                EB_FREE_ARRAY(obj->md_blk_arr_nsq[coded_leaf_index].palette_info.color_idx_map);
    }
    EB_FREE_ARRAY(obj->ref_best_ref_sq_table);
    EB_FREE_ARRAY(obj->ref_best_cost_sq_table);
    for (CandClass cand_class_it = CAND_CLASS_0; cand_class_it < CAND_CLASS_TOTAL; cand_class_it++)
        EB_FREE_ARRAY(obj->cand_buff_indices[cand_class_it]);
    EB_FREE_ARRAY(obj->best_candidate_index_array);

    EB_FREE_ARRAY(obj->above_txfm_context);
    EB_FREE_ARRAY(obj->left_txfm_context);
#if FTR_BYPASS_ENCDEC
    for (uint32_t coded_leaf_index = 0; coded_leaf_index < block_max_count_sb; ++coded_leaf_index) {
        EB_DELETE(obj->md_blk_arr_nsq[coded_leaf_index].coeff_tmp);
        EB_DELETE(obj->md_blk_arr_nsq[coded_leaf_index].recon_tmp);
}
#endif
#if NO_ENCDEC //SB128_TODO to upgrade
    int coded_leaf_index;
    for (coded_leaf_index = 0; coded_leaf_index < BLOCK_MAX_COUNT_SB_128; ++coded_leaf_index) {
        EB_DELETE(obj->md_blk_arr_nsq[coded_leaf_index].recon_tmp);
        EB_DELETE(obj->md_blk_arr_nsq[coded_leaf_index].coeff_tmp);
    }
#endif
    EB_DELETE_PTR_ARRAY(obj->candidate_buffer_ptr_array, obj->max_nics_uv);
    EB_FREE_ARRAY(obj->candidate_buffer_tx_depth_1->candidate_ptr);
    EB_DELETE(obj->candidate_buffer_tx_depth_1);
    EB_FREE_ARRAY(obj->candidate_buffer_tx_depth_2->candidate_ptr);
    EB_DELETE(obj->candidate_buffer_tx_depth_2);
    EB_DELETE(obj->trans_quant_buffers_ptr);
    EB_FREE_ALIGNED_ARRAY(obj->cfl_temp_luma_recon16bit);
    EB_FREE_ALIGNED_ARRAY(obj->cfl_temp_luma_recon);
    if (obj->is_md_rate_estimation_ptr_owner)
        EB_FREE_ARRAY(obj->md_rate_estimation_ptr);
    EB_FREE_ARRAY(obj->fast_candidate_array);
    EB_FREE_ARRAY(obj->fast_candidate_ptr_array);
    EB_FREE_ARRAY(obj->fast_cost_array);
    EB_FREE_ARRAY(obj->full_cost_array);
#if !CLN_MOVE_SKIP_MODE_CHECK
    EB_FREE_ARRAY(obj->full_cost_skip_ptr);
    EB_FREE_ARRAY(obj->full_cost_merge_ptr);
#endif
    if (obj->md_local_blk_unit) {
        EB_FREE_ARRAY(obj->md_local_blk_unit[0].neigh_left_recon_16bit[0]);
        EB_FREE_ARRAY(obj->md_local_blk_unit[0].neigh_top_recon_16bit[0]);
        EB_FREE_ARRAY(obj->md_local_blk_unit[0].neigh_left_recon[0]);
        EB_FREE_ARRAY(obj->md_local_blk_unit[0].neigh_top_recon[0]);
    }
    if (obj->md_blk_arr_nsq) {
        EB_FREE_ARRAY(obj->md_blk_arr_nsq[0].av1xd);
    }
    EB_FREE_ARRAY(obj->avail_blk_flag);
    EB_FREE_ARRAY(obj->tested_blk_flag);
    EB_FREE_ARRAY(obj->do_not_process_blk);
    EB_FREE_ARRAY(obj->md_local_blk_unit);
    EB_FREE_ARRAY(obj->md_blk_arr_nsq);
    EB_FREE_ARRAY(obj->md_ep_pipe_sb);
    EB_FREE_ARRAY(obj->mdc_sb_array);
    for (uint32_t txt_itr = 0; txt_itr < TX_TYPES; ++txt_itr) {
        EB_DELETE(obj->recon_coeff_ptr[txt_itr]);
        EB_DELETE(obj->recon_ptr[txt_itr]);
#if REFCTR_ADD_QUANT_COEFF_BUFF_MD
        EB_DELETE(obj->quant_coeff_ptr[txt_itr]);
#endif
    }
#if OPT_INTER_PRED
    EB_DELETE(obj->scratch_prediction_ptr);
#else
    EB_DELETE(obj->cfl_temp_prediction_ptr);
#endif
#if !REFCTR_ADD_QUANT_COEFF_BUFF_MD
    EB_DELETE(obj->residual_quant_coeff_ptr);
#endif
    EB_DELETE(obj->temp_residual_ptr);
    EB_DELETE(obj->temp_recon_ptr);
}


uint8_t  get_nic_scaling_level(PdPass pd_pass, EbEncMode enc_mode ,uint8_t temporal_layer_index );

#if TUNE_MDS0
/*
* return the max canidate count for MDS0
  Used by candidate injection and memory allocation
*/
uint16_t  get_max_can_count(EbEncMode enc_mode ) {
    uint16_t  max_can_count ;


    if (enc_mode <= ENC_M0)
        max_can_count = 900;
    else if (enc_mode <= ENC_M1)
        max_can_count =  720;
    else if (enc_mode <= ENC_M2)
        max_can_count = 576;
    else if (enc_mode <= ENC_M3)
        max_can_count = 461;
    else if (enc_mode <= ENC_M4)
        max_can_count = 369;
    else if (enc_mode <= ENC_M5)
        max_can_count = 295;
    else if (enc_mode <= ENC_M6)
        max_can_count = 236;
    else if (enc_mode <= ENC_M7)
        max_can_count = 190;
    else if (enc_mode <= ENC_M8)
        max_can_count = 120;
    else if (enc_mode <= ENC_M9)
        max_can_count = 120;
    else if (enc_mode <= ENC_M10)
        max_can_count = 100;
    else
        max_can_count = 80;

   return max_can_count ;
}

#endif
/******************************************************
 * Mode Decision Context Constructor
 ******************************************************/
EbErrorType mode_decision_context_ctor(ModeDecisionContext *context_ptr, EbColorFormat color_format,
                                       uint8_t sb_size,
                                        uint8_t enc_mode,
#if CLN_GEOM
                                        uint16_t max_block_cnt,
#endif
                                       EbFifo *mode_decision_configuration_input_fifo_ptr,
                                       EbFifo *mode_decision_output_fifo_ptr,
                                       uint8_t enable_hbd_mode_decision, uint8_t cfg_palette) {
    uint32_t buffer_index;
    uint32_t cand_index;

#if CLN_GEOM
    context_ptr->init_max_block_cnt = max_block_cnt;
    uint32_t block_max_count_sb = max_block_cnt;
#else
    uint32_t block_max_count_sb = (sb_size == MAX_SB_SIZE) ? BLOCK_MAX_COUNT_SB_128
                                                           : BLOCK_MAX_COUNT_SB_64;
#endif

    context_ptr->sb_size        = sb_size;
    (void)color_format;

    context_ptr->dctor             = mode_decision_context_dctor;
    context_ptr->hbd_mode_decision = enable_hbd_mode_decision;

    // Input/Output System Resource Manager FIFOs
    context_ptr->mode_decision_configuration_input_fifo_ptr =
        mode_decision_configuration_input_fifo_ptr;
    context_ptr->mode_decision_output_fifo_ptr = mode_decision_output_fifo_ptr;

    // Maximum number of candidates MD can support
    // determine MAX_NICS for a given preset
    uint32_t max_nics = 0;

    // get max number of NICS
    for (uint8_t pic_type = 0; pic_type < NICS_PIC_TYPE; pic_type++) {
        uint32_t nics = 0;
        for (CandClass cand_class_it = CAND_CLASS_0; cand_class_it < CAND_CLASS_TOTAL; cand_class_it++) {
            nics += MD_STAGE_NICS[pic_type][cand_class_it];
        }
        max_nics = MAX(max_nics ,nics);
    }


    // get the min scalling level ( smallest scalling level is the most concervative
    uint8_t min_nic_scaling_level = NICS_SCALING_LEVELS - 1 ;
    for (PdPass pd_id = 0; pd_id < PD_PASS_TOTAL; pd_id++) {

        for (uint8_t temporal_layer_index = 0; temporal_layer_index < MAX_TEMPORAL_LAYERS; temporal_layer_index++) {
            min_nic_scaling_level = MIN(min_nic_scaling_level, get_nic_scaling_level(pd_id,  enc_mode,  temporal_layer_index));
        }
    }
    // scale max_nics
    uint8_t stage1_scaling_num = MD_STAGE_NICS_SCAL_NUM[min_nic_scaling_level][MD_STAGE_1];
    max_nics   = MAX(2,
            DIVIDE_AND_ROUND(max_nics * stage1_scaling_num, MD_STAGE_NICS_SCAL_DENUM));

    max_nics += CAND_CLASS_TOTAL ; //need one extra temp buffer for each fast loop call
    context_ptr->max_nics = max_nics;
    context_ptr->max_nics_uv = max_nics + 84; // needed for independant chroma search
    // Cfl scratch memory
    if (context_ptr->hbd_mode_decision > EB_8_BIT_MD)
        EB_MALLOC_ALIGNED(context_ptr->cfl_temp_luma_recon16bit,
                          sizeof(uint16_t) * sb_size * sb_size);
    if (context_ptr->hbd_mode_decision != EB_10_BIT_MD)
        EB_MALLOC_ALIGNED(context_ptr->cfl_temp_luma_recon, sizeof(uint8_t) * sb_size * sb_size);
    // MD rate Estimation tables
    EB_MALLOC_ARRAY(context_ptr->md_rate_estimation_ptr, 1);
    context_ptr->is_md_rate_estimation_ptr_owner = EB_TRUE;

    EB_MALLOC_ARRAY(context_ptr->md_local_blk_unit, block_max_count_sb);
    EB_MALLOC_ARRAY(context_ptr->md_blk_arr_nsq, block_max_count_sb);
    EB_MALLOC_ARRAY(context_ptr->md_ep_pipe_sb, block_max_count_sb);
    // Fast Candidate Array
#if TUNE_MDS0
    uint16_t max_can_count = get_max_can_count (enc_mode) + 84 ;//chroma search;
    EB_MALLOC_ARRAY(context_ptr->fast_candidate_array, max_can_count);

    EB_MALLOC_ARRAY(context_ptr->fast_candidate_ptr_array, max_can_count);

    for (cand_index = 0; cand_index < max_can_count; ++cand_index) {
        context_ptr->fast_candidate_ptr_array[cand_index] =
            &context_ptr->fast_candidate_array[cand_index];
        context_ptr->fast_candidate_ptr_array[cand_index]->palette_info = NULL;
    }
#else
    EB_MALLOC_ARRAY(context_ptr->fast_candidate_array, MODE_DECISION_CANDIDATE_MAX_COUNT);

    EB_MALLOC_ARRAY(context_ptr->fast_candidate_ptr_array, MODE_DECISION_CANDIDATE_MAX_COUNT);

    for (cand_index = 0; cand_index < MODE_DECISION_CANDIDATE_MAX_COUNT; ++cand_index) {
        context_ptr->fast_candidate_ptr_array[cand_index] =
            &context_ptr->fast_candidate_array[cand_index];
        context_ptr->fast_candidate_ptr_array[cand_index]->palette_info = NULL;
    }
#endif

    for (int cd = 0; cd < MAX_PAL_CAND; cd++)
        if (cfg_palette)
            EB_MALLOC_ARRAY(context_ptr->palette_cand_array[cd].color_idx_map, MAX_PALETTE_SQUARE);
        else
            context_ptr->palette_cand_array[cd].color_idx_map = NULL;
    // Transform and Quantization Buffers
    EB_NEW(context_ptr->trans_quant_buffers_ptr, svt_trans_quant_buffers_ctor, sb_size);

    // Cost Arrays
    EB_MALLOC_ARRAY(context_ptr->fast_cost_array,  context_ptr->max_nics_uv);
    EB_MALLOC_ARRAY(context_ptr->full_cost_array,  context_ptr->max_nics_uv);
#if !CLN_MOVE_SKIP_MODE_CHECK
    EB_MALLOC_ARRAY(context_ptr->full_cost_skip_ptr,  context_ptr->max_nics_uv);
    EB_MALLOC_ARRAY(context_ptr->full_cost_merge_ptr,  context_ptr->max_nics_uv);
#endif
    // Candidate Buffers
    EB_NEW(context_ptr->candidate_buffer_tx_depth_1,
           mode_decision_scratch_candidate_buffer_ctor,
           sb_size,
           context_ptr->hbd_mode_decision ? EB_10BIT : EB_8BIT);

    EB_ALLOC_PTR_ARRAY(context_ptr->candidate_buffer_tx_depth_1->candidate_ptr, 1);
    EB_NEW(context_ptr->candidate_buffer_tx_depth_2,
           mode_decision_scratch_candidate_buffer_ctor,
           sb_size,
           context_ptr->hbd_mode_decision ? EB_10BIT : EB_8BIT);

    EB_ALLOC_PTR_ARRAY(context_ptr->candidate_buffer_tx_depth_2->candidate_ptr, 1);
    context_ptr->md_local_blk_unit[0].neigh_left_recon[0]       = NULL;
    context_ptr->md_local_blk_unit[0].neigh_top_recon[0]        = NULL;
    context_ptr->md_local_blk_unit[0].neigh_left_recon_16bit[0] = NULL;
    context_ptr->md_local_blk_unit[0].neigh_top_recon_16bit[0]  = NULL;
    uint16_t sz                                                 = sizeof(uint16_t);
    if (context_ptr->hbd_mode_decision > EB_8_BIT_MD) {
        EB_MALLOC_ARRAY(context_ptr->md_local_blk_unit[0].neigh_left_recon_16bit[0],
                        block_max_count_sb * sb_size * 3 * sz);
        EB_MALLOC_ARRAY(context_ptr->md_local_blk_unit[0].neigh_top_recon_16bit[0],
                        block_max_count_sb * sb_size * 3 * sz);
    }
    if (context_ptr->hbd_mode_decision != EB_10_BIT_MD) {
        EB_MALLOC_ARRAY(context_ptr->md_local_blk_unit[0].neigh_left_recon[0],
                        block_max_count_sb * sb_size * 3);
        EB_MALLOC_ARRAY(context_ptr->md_local_blk_unit[0].neigh_top_recon[0],
                        block_max_count_sb * sb_size * 3);
    }
    uint32_t coded_leaf_index;
    for (coded_leaf_index = 0; coded_leaf_index < block_max_count_sb; ++coded_leaf_index) {
        for (int i = 0; i < 3; i++) {
            size_t offset = (coded_leaf_index * sb_size * 3 + i * sb_size) * sz;
            context_ptr->md_local_blk_unit[coded_leaf_index].neigh_left_recon_16bit[i] =
                context_ptr->md_local_blk_unit[0].neigh_left_recon_16bit[0] + offset;
            context_ptr->md_local_blk_unit[coded_leaf_index].neigh_top_recon_16bit[i] =
                context_ptr->md_local_blk_unit[0].neigh_top_recon_16bit[0] + offset;
        }
        for (int i = 0; i < 3; i++) {
            size_t offset = coded_leaf_index * sb_size * 3 + i * sb_size;
            context_ptr->md_local_blk_unit[coded_leaf_index].neigh_left_recon[i] =
                context_ptr->md_local_blk_unit[0].neigh_left_recon[0] + offset;
            context_ptr->md_local_blk_unit[coded_leaf_index].neigh_top_recon[i] =
                context_ptr->md_local_blk_unit[0].neigh_top_recon[0] + offset;
        }
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
        if (svt_av1_allow_palette(cfg_palette, blk_geom->bsize))
            EB_MALLOC_ARRAY(
                context_ptr->md_blk_arr_nsq[coded_leaf_index].palette_info.color_idx_map,
                MAX_PALETTE_SQUARE);
        else
            context_ptr->md_blk_arr_nsq[coded_leaf_index].palette_info.color_idx_map = NULL;
#if FTR_BYPASS_ENCDEC
        {
            EbPictureBufferDescInitData init_data;

            init_data.buffer_enable_mask    = PICTURE_BUFFER_DESC_FULL_MASK;
            init_data.max_width             = blk_geom->bwidth;
            init_data.max_height            = blk_geom->bheight;
            init_data.bit_depth             = EB_32BIT;
            init_data.color_format          = (blk_geom->bwidth > 4 && blk_geom->bheight > 4) ? EB_YUV420 : EB_YUV444; // PW - must have at least 4x4 for chroma coeffs
            init_data.left_padding          = 0;
            init_data.right_padding         = 0;
            init_data.top_padding           = 0;
            init_data.bot_padding           = 0;
            init_data.split_mode            = EB_FALSE;

            EB_NEW(context_ptr->md_blk_arr_nsq[coded_leaf_index].coeff_tmp,
                svt_picture_buffer_desc_ctor,
                (EbPtr)&init_data);

            init_data.buffer_enable_mask    = PICTURE_BUFFER_DESC_FULL_MASK;
            init_data.max_width             = blk_geom->bwidth;
            init_data.max_height            = blk_geom->bheight;
            init_data.bit_depth             = context_ptr->hbd_mode_decision ? EB_10BIT : EB_8BIT;;
            init_data.color_format          = (blk_geom->bwidth > 4 && blk_geom->bheight > 4) ? EB_YUV420 : EB_YUV444;
            init_data.left_padding          = 0;
            init_data.right_padding         = 0;
            init_data.top_padding           = 0;
            init_data.bot_padding           = 0;
            init_data.split_mode            = EB_FALSE;

            EB_NEW(context_ptr->md_blk_arr_nsq[coded_leaf_index].recon_tmp,
                svt_picture_buffer_desc_ctor,
                (EbPtr)&init_data);
    }
#endif
#if NO_ENCDEC //SB128_TODO to upgrade
        {
            EbPictureBufferDescInitData init_data;

            init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
            init_data.max_width          = SB_STRIDE_Y;
            init_data.max_height         = SB_STRIDE_Y;
            init_data.bit_depth          = EB_32BIT;
            init_data.color_format       = EB_YUV420;
            init_data.left_padding       = 0;
            init_data.right_padding      = 0;
            init_data.top_padding        = 0;
            init_data.bot_padding        = 0;
            init_data.split_mode         = EB_FALSE;

            EB_NEW(context_ptr->md_blk_arr_nsq[coded_leaf_index].coeff_tmp,
                   svt_picture_buffer_desc_ctor,
                   (EbPtr)&init_data);

            init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
            init_data.max_width          = SB_STRIDE_Y;
            init_data.max_height         = SB_STRIDE_Y;
            init_data.bit_depth          = EB_8BIT;
            init_data.color_format       = EB_YUV420;
            init_data.left_padding       = 0;
            init_data.right_padding      = 0;
            init_data.top_padding        = 0;
            init_data.bot_padding        = 0;
            init_data.split_mode         = EB_FALSE;

            EB_NEW(context_ptr->md_blk_arr_nsq[coded_leaf_index].recon_tmp,
                   svt_picture_buffer_desc_ctor,
                   (EbPtr)&init_data);
        }
#endif
    }
    for (CandClass cand_class_it = CAND_CLASS_0; cand_class_it < CAND_CLASS_TOTAL; cand_class_it++)
        EB_MALLOC_ARRAY(context_ptr->cand_buff_indices[cand_class_it], context_ptr->max_nics_uv);

    EB_MALLOC_ARRAY(context_ptr->best_candidate_index_array, context_ptr->max_nics_uv);

    EB_MALLOC_ARRAY(context_ptr->ref_best_cost_sq_table, MAX_REF_TYPE_CAND);
    EB_MALLOC_ARRAY(context_ptr->ref_best_ref_sq_table, MAX_REF_TYPE_CAND);
    EB_MALLOC_ARRAY(context_ptr->above_txfm_context, (sb_size >> MI_SIZE_LOG2));
    EB_MALLOC_ARRAY(context_ptr->left_txfm_context, (sb_size >> MI_SIZE_LOG2));
    EbPictureBufferDescInitData thirty_two_width_picture_buffer_desc_init_data;
    EbPictureBufferDescInitData picture_buffer_desc_init_data;

    picture_buffer_desc_init_data.max_width  = sb_size;
    picture_buffer_desc_init_data.max_height = sb_size;
    picture_buffer_desc_init_data.bit_depth  = context_ptr->hbd_mode_decision ? EB_10BIT : EB_8BIT;
    picture_buffer_desc_init_data.color_format       = EB_YUV420;
    picture_buffer_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    picture_buffer_desc_init_data.left_padding       = 0;
    picture_buffer_desc_init_data.right_padding      = 0;
    picture_buffer_desc_init_data.top_padding        = 0;
    picture_buffer_desc_init_data.bot_padding        = 0;
    picture_buffer_desc_init_data.split_mode         = EB_FALSE;

    thirty_two_width_picture_buffer_desc_init_data.max_width    = sb_size;
    thirty_two_width_picture_buffer_desc_init_data.max_height   = sb_size;
    thirty_two_width_picture_buffer_desc_init_data.bit_depth    = EB_32BIT;
    thirty_two_width_picture_buffer_desc_init_data.color_format = EB_YUV420;
    thirty_two_width_picture_buffer_desc_init_data.buffer_enable_mask =
        PICTURE_BUFFER_DESC_FULL_MASK;
    thirty_two_width_picture_buffer_desc_init_data.left_padding  = 0;
    thirty_two_width_picture_buffer_desc_init_data.right_padding = 0;
    thirty_two_width_picture_buffer_desc_init_data.top_padding   = 0;
    thirty_two_width_picture_buffer_desc_init_data.bot_padding   = 0;
    thirty_two_width_picture_buffer_desc_init_data.split_mode    = EB_FALSE;

    for (uint32_t txt_itr = 0; txt_itr < TX_TYPES; ++txt_itr) {
        EB_NEW(context_ptr->recon_coeff_ptr[txt_itr],
               svt_picture_buffer_desc_ctor,
               (EbPtr)&thirty_two_width_picture_buffer_desc_init_data);
        EB_NEW(context_ptr->recon_ptr[txt_itr],
               svt_picture_buffer_desc_ctor,
               (EbPtr)&picture_buffer_desc_init_data);
#if REFCTR_ADD_QUANT_COEFF_BUFF_MD
        EB_NEW(context_ptr->quant_coeff_ptr[txt_itr],
            svt_picture_buffer_desc_ctor,
            (EbPtr)&thirty_two_width_picture_buffer_desc_init_data);
#endif
    }
#if !REFCTR_ADD_QUANT_COEFF_BUFF_MD
    EB_NEW(context_ptr->residual_quant_coeff_ptr,
           svt_picture_buffer_desc_ctor,
           (EbPtr)&thirty_two_width_picture_buffer_desc_init_data);
#endif
#if OPT_INTER_PRED
    EB_NEW(context_ptr->scratch_prediction_ptr,
        svt_picture_buffer_desc_ctor,
        (EbPtr)&picture_buffer_desc_init_data);
#else
    EB_NEW(context_ptr->cfl_temp_prediction_ptr,
           svt_picture_buffer_desc_ctor,
           (EbPtr)&picture_buffer_desc_init_data);
#endif
    EbPictureBufferDescInitData double_width_picture_buffer_desc_init_data;
    double_width_picture_buffer_desc_init_data.max_width          = sb_size;
    double_width_picture_buffer_desc_init_data.max_height         = sb_size;
    double_width_picture_buffer_desc_init_data.bit_depth          = EB_16BIT;
    double_width_picture_buffer_desc_init_data.color_format       = EB_YUV420;
    double_width_picture_buffer_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    double_width_picture_buffer_desc_init_data.left_padding       = 0;
    double_width_picture_buffer_desc_init_data.right_padding      = 0;
    double_width_picture_buffer_desc_init_data.top_padding        = 0;
    double_width_picture_buffer_desc_init_data.bot_padding        = 0;
    double_width_picture_buffer_desc_init_data.split_mode         = EB_FALSE;

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

    for (buffer_index = 0; buffer_index <  context_ptr->max_nics; ++buffer_index) {
#if CLN_MOVE_SKIP_MODE_CHECK
        EB_NEW(context_ptr->candidate_buffer_ptr_array[buffer_index],
            mode_decision_candidate_buffer_ctor,
            context_ptr->hbd_mode_decision ? EB_10BIT : EB_8BIT,
            sb_size,
            PICTURE_BUFFER_DESC_FULL_MASK,
            context_ptr->temp_residual_ptr,
            context_ptr->temp_recon_ptr,
            &(context_ptr->fast_cost_array[buffer_index]),
            &(context_ptr->full_cost_array[buffer_index]));
#else
        EB_NEW(context_ptr->candidate_buffer_ptr_array[buffer_index],
               mode_decision_candidate_buffer_ctor,
               context_ptr->hbd_mode_decision ? EB_10BIT : EB_8BIT,
               sb_size,
               PICTURE_BUFFER_DESC_FULL_MASK,
               context_ptr->temp_residual_ptr,
               context_ptr->temp_recon_ptr,
               &(context_ptr->fast_cost_array[buffer_index]),
               &(context_ptr->full_cost_array[buffer_index]),
               &(context_ptr->full_cost_skip_ptr[buffer_index]),
               &(context_ptr->full_cost_merge_ptr[buffer_index]));
#endif
    }

    for (buffer_index = max_nics; buffer_index <  context_ptr->max_nics_uv; ++buffer_index) {
#if CLN_MOVE_SKIP_MODE_CHECK
        EB_NEW(context_ptr->candidate_buffer_ptr_array[buffer_index],
            mode_decision_candidate_buffer_ctor,
            context_ptr->hbd_mode_decision ? EB_10BIT : EB_8BIT,
            sb_size,
            PICTURE_BUFFER_DESC_CHROMA_MASK,
            context_ptr->temp_residual_ptr,
            context_ptr->temp_recon_ptr,
            &(context_ptr->fast_cost_array[buffer_index]),
            &(context_ptr->full_cost_array[buffer_index]));
#else
        EB_NEW(context_ptr->candidate_buffer_ptr_array[buffer_index],
               mode_decision_candidate_buffer_ctor,
               context_ptr->hbd_mode_decision ? EB_10BIT : EB_8BIT,
               sb_size,
               PICTURE_BUFFER_DESC_CHROMA_MASK,
               context_ptr->temp_residual_ptr,
               context_ptr->temp_recon_ptr,
               &(context_ptr->fast_cost_array[buffer_index]),
               &(context_ptr->full_cost_array[buffer_index]),
               &(context_ptr->full_cost_skip_ptr[buffer_index]),
               &(context_ptr->full_cost_merge_ptr[buffer_index]));
#endif
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
        neighbor_array_unit_reset(pcs_ptr->md_ref_frame_type_neighbor_array[depth][tile_idx]);

        neighbor_array_unit_reset32(pcs_ptr->md_interpolation_type_neighbor_array[depth][tile_idx]);
    }

    return;
}

// Set the lambda for each sb.
// When lambda tuning is on (blk_lambda_tuning), lambda of each block is set separately (full_lambda_md/fast_lambda_md)
// later in set_tuned_blk_lambda
// Testing showed that updating SAD lambda based on frame info was not helpful; therefore, the SAD lambda generation is not changed.
int compute_rdmult_sse(PictureControlSet *pcs_ptr, uint8_t q_index, uint8_t bit_depth);

void av1_lambda_assign_md(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr) {
    context_ptr->full_lambda_md[0] = (uint32_t)compute_rdmult_sse(
        pcs_ptr, context_ptr->qp_index, 8);
    context_ptr->fast_lambda_md[0] = av1_lambda_mode_decision8_bit_sad[context_ptr->qp_index];

    context_ptr->full_lambda_md[1] = (uint32_t)compute_rdmult_sse(
        pcs_ptr, context_ptr->qp_index, 10);
    context_ptr->fast_lambda_md[1] = av1lambda_mode_decision10_bit_sad[context_ptr->qp_index];

    context_ptr->full_lambda_md[1] *= 16;
    context_ptr->fast_lambda_md[1] *= 4;
    context_ptr->full_sb_lambda_md[0] = context_ptr->full_lambda_md[0];
    context_ptr->full_sb_lambda_md[1] = context_ptr->full_lambda_md[1];
}

void av1_lambda_assign(PictureControlSet *pcs_ptr, uint32_t *fast_lambda, uint32_t *full_lambda,
                       uint8_t bit_depth, uint16_t qp_index, EbBool multiply_lambda) {
    if (bit_depth == 8) {
        *full_lambda = (uint32_t)compute_rdmult_sse(pcs_ptr, (uint8_t)qp_index, bit_depth);
        *fast_lambda = av1_lambda_mode_decision8_bit_sad[qp_index];
    } else if (bit_depth == 10) {
        *full_lambda = (uint32_t)compute_rdmult_sse(pcs_ptr, (uint8_t)qp_index, bit_depth);
        *fast_lambda = av1lambda_mode_decision10_bit_sad[qp_index];
        if (multiply_lambda) {
            *full_lambda *= 16;
            *fast_lambda *= 4;
        }
    } else if (bit_depth == 12) {
        *full_lambda = (uint32_t)compute_rdmult_sse(pcs_ptr, (uint8_t)qp_index, bit_depth);
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
    FrameHeader *frm_hdr           = &pcs_ptr->parent_pcs_ptr->frm_hdr;
    context_ptr->hbd_mode_decision = pcs_ptr->hbd_mode_decision;
    // QP
    context_ptr->qp_index = (uint8_t)frm_hdr->quantization_params.base_q_idx;
    av1_lambda_assign_md(pcs_ptr, context_ptr);
    // Reset MD rate Estimation table to initial values by copying from md_rate_estimation_array
    if (context_ptr->is_md_rate_estimation_ptr_owner) {
        context_ptr->is_md_rate_estimation_ptr_owner = EB_FALSE;
        EB_FREE_ARRAY(context_ptr->md_rate_estimation_ptr);
    }
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
    return;
}

/******************************************************
 * Mode Decision Configure SB
 ******************************************************/
void mode_decision_configure_sb(ModeDecisionContext *context_ptr, PictureControlSet *pcs_ptr,
                                uint8_t sb_qp) {
    /* Note(CHKN) : when Qp modulation varies QP on a sub-SB(CU) basis,  Lamda has to change based on Cu->QP , and then this code has to move inside the CU loop in MD */

    // Lambda Assignement
    context_ptr->qp_index = pcs_ptr->parent_pcs_ptr->frm_hdr.delta_q_params.delta_q_present
        ? sb_qp
        : (uint8_t)pcs_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx;

    av1_lambda_assign_md(pcs_ptr, context_ptr);

    return;
}
