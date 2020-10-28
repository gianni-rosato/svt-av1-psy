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

int svt_av1_allow_palette(int allow_palette, BlockSize sb_type);
static void mode_decision_context_dctor(EbPtr p) {
    ModeDecisionContext *obj = (ModeDecisionContext *)p;
    uint32_t block_max_count_sb = (obj->sb_size == MAX_SB_SIZE) ? BLOCK_MAX_COUNT_SB_128 :
                                                                  BLOCK_MAX_COUNT_SB_64;
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
    EB_FREE_ARRAY(obj->above_txfm_context);
    EB_FREE_ARRAY(obj->left_txfm_context);
#if NO_ENCDEC //SB128_TODO to upgrade
    int coded_leaf_index;
    for (coded_leaf_index = 0; coded_leaf_index < BLOCK_MAX_COUNT_SB_128; ++coded_leaf_index) {
        EB_DELETE(obj->md_blk_arr_nsq[coded_leaf_index].recon_tmp);
        EB_DELETE(obj->md_blk_arr_nsq[coded_leaf_index].coeff_tmp);
    }
#endif
    EB_DELETE_PTR_ARRAY(obj->candidate_buffer_ptr_array, MAX_NFL_BUFF);
    EB_FREE_ARRAY(obj->candidate_buffer_tx_depth_1->candidate_ptr);
    EB_DELETE(obj->candidate_buffer_tx_depth_1);
    EB_FREE_ARRAY(obj->candidate_buffer_tx_depth_2->candidate_ptr);
    EB_DELETE(obj->candidate_buffer_tx_depth_2);
    EB_DELETE(obj->trans_quant_buffers_ptr);
    EB_FREE_ALIGNED_ARRAY(obj->cfl_temp_luma_recon16bit);
    EB_FREE_ALIGNED_ARRAY(obj->cfl_temp_luma_recon);
    if (obj->is_md_rate_estimation_ptr_owner) EB_FREE_ARRAY(obj->md_rate_estimation_ptr);
    EB_FREE_ARRAY(obj->fast_candidate_array);
    EB_FREE_ARRAY(obj->fast_candidate_ptr_array);
    EB_FREE_ARRAY(obj->fast_cost_array);
    EB_FREE_ARRAY(obj->full_cost_array);
    EB_FREE_ARRAY(obj->full_cost_skip_ptr);
    EB_FREE_ARRAY(obj->full_cost_merge_ptr);
    if (obj->md_local_blk_unit) {
        EB_FREE_ARRAY(obj->md_local_blk_unit[0].neigh_left_recon_16bit[0]);
        EB_FREE_ARRAY(obj->md_local_blk_unit[0].neigh_top_recon_16bit[0]);
        EB_FREE_ARRAY(obj->md_local_blk_unit[0].neigh_left_recon[0]);
        EB_FREE_ARRAY(obj->md_local_blk_unit[0].neigh_top_recon[0]);
    }
    if (obj->md_blk_arr_nsq) {
        EB_FREE_ARRAY(obj->md_blk_arr_nsq[0].av1xd);
    }
    EB_FREE_ARRAY(obj->md_local_blk_unit);
    EB_FREE_ARRAY(obj->md_blk_arr_nsq);
    EB_FREE_ARRAY(obj->md_ep_pipe_sb);
    EB_FREE_ARRAY(obj->mdc_sb_array);
    for (uint32_t txt_itr = 0; txt_itr < TX_TYPES; ++txt_itr) {
        EB_DELETE(obj->recon_coeff_ptr[txt_itr]);
        EB_DELETE(obj->recon_ptr[txt_itr]);
    }
    EB_DELETE(obj->prediction_ptr_temp);
    EB_DELETE(obj->cfl_temp_prediction_ptr);
    EB_DELETE(obj->residual_quant_coeff_ptr);

    EB_DELETE(obj->temp_residual_ptr);
    EB_DELETE(obj->temp_recon_ptr);
}

/******************************************************
 * Mode Decision Context Constructor
 ******************************************************/
EbErrorType mode_decision_context_ctor(ModeDecisionContext *context_ptr, EbColorFormat color_format,
                                       uint8_t sb_size,
                                       EbFifo *mode_decision_configuration_input_fifo_ptr,
                                       EbFifo *mode_decision_output_fifo_ptr,
                                       uint8_t enable_hbd_mode_decision, uint8_t cfg_palette) {
    uint32_t buffer_index;
    uint32_t cand_index;
    uint32_t block_max_count_sb = (sb_size == MAX_SB_SIZE) ? BLOCK_MAX_COUNT_SB_128 :
                                                             BLOCK_MAX_COUNT_SB_64;
    context_ptr->sb_size = sb_size;
    (void)color_format;

    context_ptr->dctor             = mode_decision_context_dctor;
    context_ptr->hbd_mode_decision = enable_hbd_mode_decision;

    // Input/Output System Resource Manager FIFOs
    context_ptr->mode_decision_configuration_input_fifo_ptr =
        mode_decision_configuration_input_fifo_ptr;
    context_ptr->mode_decision_output_fifo_ptr = mode_decision_output_fifo_ptr;

    // Cfl scratch memory
    if (context_ptr->hbd_mode_decision > EB_8_BIT_MD)
        EB_MALLOC_ALIGNED(context_ptr->cfl_temp_luma_recon16bit, sizeof(uint16_t) * sb_size * sb_size);
    if (context_ptr->hbd_mode_decision != EB_10_BIT_MD)
        EB_MALLOC_ALIGNED(context_ptr->cfl_temp_luma_recon, sizeof(uint8_t) * sb_size * sb_size);
    // MD rate Estimation tables
    EB_MALLOC_ARRAY(context_ptr->md_rate_estimation_ptr, 1);
    context_ptr->is_md_rate_estimation_ptr_owner = EB_TRUE;

    EB_MALLOC_ARRAY(context_ptr->md_local_blk_unit, block_max_count_sb);
    EB_MALLOC_ARRAY(context_ptr->md_blk_arr_nsq, block_max_count_sb);
    EB_MALLOC_ARRAY(context_ptr->md_ep_pipe_sb, block_max_count_sb);
    // Fast Candidate Array
    EB_MALLOC_ARRAY(context_ptr->fast_candidate_array, MODE_DECISION_CANDIDATE_MAX_COUNT);

    EB_MALLOC_ARRAY(context_ptr->fast_candidate_ptr_array, MODE_DECISION_CANDIDATE_MAX_COUNT);

    for (cand_index = 0; cand_index < MODE_DECISION_CANDIDATE_MAX_COUNT; ++cand_index) {
        context_ptr->fast_candidate_ptr_array[cand_index] =
            &context_ptr->fast_candidate_array[cand_index];
        context_ptr->fast_candidate_ptr_array[cand_index]->md_rate_estimation_ptr =
            context_ptr->md_rate_estimation_ptr;
            context_ptr->fast_candidate_ptr_array[cand_index]->palette_info = NULL;
    }
    for (int cd = 0; cd < MAX_PAL_CAND; cd++)
        if (cfg_palette)
            EB_MALLOC_ARRAY(context_ptr->palette_cand_array[cd].color_idx_map, MAX_PALETTE_SQUARE);
        else
            context_ptr->palette_cand_array[cd].color_idx_map = NULL;
    // Transform and Quantization Buffers
    EB_NEW(context_ptr->trans_quant_buffers_ptr, svt_trans_quant_buffers_ctor, sb_size);

    // Cost Arrays
    EB_MALLOC_ARRAY(context_ptr->fast_cost_array, MAX_NFL_BUFF);
    EB_MALLOC_ARRAY(context_ptr->full_cost_array, MAX_NFL_BUFF);
    EB_MALLOC_ARRAY(context_ptr->full_cost_skip_ptr, MAX_NFL_BUFF);
    EB_MALLOC_ARRAY(context_ptr->full_cost_merge_ptr, MAX_NFL_BUFF);
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
    context_ptr->md_blk_arr_nsq[0].av1xd                     = NULL;
    EB_MALLOC_ARRAY(context_ptr->md_blk_arr_nsq[0].av1xd, block_max_count_sb);

    EB_MALLOC_ARRAY(context_ptr->mdc_sb_array, 1);
    for (coded_leaf_index = 0; coded_leaf_index < block_max_count_sb; ++coded_leaf_index) {
        context_ptr->md_blk_arr_nsq[coded_leaf_index].av1xd =
            context_ptr->md_blk_arr_nsq[0].av1xd + coded_leaf_index;
        context_ptr->md_blk_arr_nsq[coded_leaf_index].segment_id = 0;
        const BlockGeom *blk_geom = get_blk_geom_mds(coded_leaf_index);
        if (svt_av1_allow_palette(cfg_palette, blk_geom->bsize))
            EB_MALLOC_ARRAY(
                context_ptr->md_blk_arr_nsq[coded_leaf_index].palette_info.color_idx_map,
                MAX_PALETTE_SQUARE);
        else
            context_ptr->md_blk_arr_nsq[coded_leaf_index].palette_info.color_idx_map = NULL;
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
    EB_MALLOC_ARRAY(context_ptr->ref_best_cost_sq_table, MAX_REF_TYPE_CAND);
    EB_MALLOC_ARRAY(context_ptr->ref_best_ref_sq_table, MAX_REF_TYPE_CAND);
    EB_MALLOC_ARRAY(context_ptr->above_txfm_context, (sb_size >> MI_SIZE_LOG2));
    EB_MALLOC_ARRAY(context_ptr->left_txfm_context, (sb_size >> MI_SIZE_LOG2));
    EbPictureBufferDescInitData thirty_two_width_picture_buffer_desc_init_data;
    EbPictureBufferDescInitData picture_buffer_desc_init_data;

    picture_buffer_desc_init_data.max_width = sb_size;
    picture_buffer_desc_init_data.max_height = sb_size;
    picture_buffer_desc_init_data.bit_depth = context_ptr->hbd_mode_decision ? EB_10BIT : EB_8BIT;
    picture_buffer_desc_init_data.color_format = EB_YUV420;
    picture_buffer_desc_init_data.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
    picture_buffer_desc_init_data.left_padding = 0;
    picture_buffer_desc_init_data.right_padding = 0;
    picture_buffer_desc_init_data.top_padding = 0;
    picture_buffer_desc_init_data.bot_padding = 0;
    picture_buffer_desc_init_data.split_mode = EB_FALSE;

    thirty_two_width_picture_buffer_desc_init_data.max_width = sb_size;
    thirty_two_width_picture_buffer_desc_init_data.max_height = sb_size;
    thirty_two_width_picture_buffer_desc_init_data.bit_depth = EB_32BIT;
    thirty_two_width_picture_buffer_desc_init_data.color_format = EB_YUV420;
    thirty_two_width_picture_buffer_desc_init_data.buffer_enable_mask =
        PICTURE_BUFFER_DESC_FULL_MASK;
    thirty_two_width_picture_buffer_desc_init_data.left_padding = 0;
    thirty_two_width_picture_buffer_desc_init_data.right_padding = 0;
    thirty_two_width_picture_buffer_desc_init_data.top_padding = 0;
    thirty_two_width_picture_buffer_desc_init_data.bot_padding = 0;
    thirty_two_width_picture_buffer_desc_init_data.split_mode = EB_FALSE;

    for (uint32_t txt_itr = 0; txt_itr < TX_TYPES; ++txt_itr) {
        EB_NEW(context_ptr->recon_coeff_ptr[txt_itr],
            svt_picture_buffer_desc_ctor,
            (EbPtr)&thirty_two_width_picture_buffer_desc_init_data);
        EB_NEW(context_ptr->recon_ptr[txt_itr],
            svt_picture_buffer_desc_ctor,
            (EbPtr)&picture_buffer_desc_init_data);
    }
    EB_NEW(context_ptr->residual_quant_coeff_ptr,
        svt_picture_buffer_desc_ctor,
        (EbPtr)&thirty_two_width_picture_buffer_desc_init_data);

    EB_NEW(context_ptr->prediction_ptr_temp,
        svt_picture_buffer_desc_ctor,
        (EbPtr)&picture_buffer_desc_init_data);

    EB_NEW(context_ptr->cfl_temp_prediction_ptr,
        svt_picture_buffer_desc_ctor,
        (EbPtr)&picture_buffer_desc_init_data);

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
    EB_ALLOC_PTR_ARRAY(context_ptr->candidate_buffer_ptr_array, MAX_NFL_BUFF);
    for (buffer_index = 0; buffer_index < MAX_NFL_BUFF_Y; ++buffer_index) {
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
    }

    for (buffer_index = MAX_NFL_BUFF_Y; buffer_index < MAX_NFL_BUFF; ++buffer_index) {
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
        neighbor_array_unit_reset(pcs_ptr->md_intra_chroma_mode_neighbor_array[depth][tile_idx]);
        neighbor_array_unit_reset(pcs_ptr->md_mv_neighbor_array[depth][tile_idx]);
        neighbor_array_unit_reset(pcs_ptr->md_skip_flag_neighbor_array[depth][tile_idx]);
        neighbor_array_unit_reset(pcs_ptr->md_mode_type_neighbor_array[depth][tile_idx]);
        neighbor_array_unit_reset(pcs_ptr->md_leaf_depth_neighbor_array[depth][tile_idx]);
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

        neighbor_array_unit_reset(pcs_ptr->md_skip_coeff_neighbor_array[depth][tile_idx]);
        neighbor_array_unit_reset(
            pcs_ptr->md_luma_dc_sign_level_coeff_neighbor_array[depth][tile_idx]);
        neighbor_array_unit_reset(
            pcs_ptr->md_tx_depth_1_luma_dc_sign_level_coeff_neighbor_array[depth][tile_idx]);
        neighbor_array_unit_reset(
            pcs_ptr->md_cb_dc_sign_level_coeff_neighbor_array[depth][tile_idx]);
        neighbor_array_unit_reset(
            pcs_ptr->md_cr_dc_sign_level_coeff_neighbor_array[depth][tile_idx]);
        neighbor_array_unit_reset(pcs_ptr->md_txfm_context_array[depth][tile_idx]);
        neighbor_array_unit_reset(pcs_ptr->md_inter_pred_dir_neighbor_array[depth][tile_idx]);
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

void av1_lambda_assign_md(PictureControlSet *pcs_ptr,
                          ModeDecisionContext   *context_ptr)
{
        context_ptr->full_lambda_md[0] = (uint32_t)compute_rdmult_sse(pcs_ptr, context_ptr->qp_index, 8);
        context_ptr->fast_lambda_md[0] = av1_lambda_mode_decision8_bit_sad[context_ptr->qp_index];

        context_ptr->full_lambda_md[1] = (uint32_t)compute_rdmult_sse(pcs_ptr, context_ptr->qp_index, 10);
        context_ptr->fast_lambda_md[1] = av1lambda_mode_decision10_bit_sad[context_ptr->qp_index];

        context_ptr->full_lambda_md[1] *= 16;
        context_ptr->fast_lambda_md[1] *= 4;
        context_ptr->full_sb_lambda_md[0] = context_ptr->full_lambda_md[0];
        context_ptr->full_sb_lambda_md[1] = context_ptr->full_lambda_md[1];
}

void av1_lambda_assign(PictureControlSet *pcs_ptr, uint32_t *fast_lambda, uint32_t *full_lambda, uint8_t bit_depth, uint16_t qp_index,
                       EbBool multiply_lambda) {
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
    FrameHeader *frm_hdr = &pcs_ptr->parent_pcs_ptr->frm_hdr;
    context_ptr->hbd_mode_decision = pcs_ptr->hbd_mode_decision;
    // QP
    context_ptr->qp_index  = (uint8_t)frm_hdr->quantization_params.base_q_idx;
    av1_lambda_assign_md(pcs_ptr, context_ptr);
    // Reset MD rate Estimation table to initial values by copying from md_rate_estimation_array
    if (context_ptr->is_md_rate_estimation_ptr_owner) {
        context_ptr->is_md_rate_estimation_ptr_owner = EB_FALSE;
        EB_FREE_ARRAY(context_ptr->md_rate_estimation_ptr);
    }
    context_ptr->md_rate_estimation_ptr = pcs_ptr->md_rate_estimation_array;
    uint32_t cand_index;
    for (cand_index = 0; cand_index < MODE_DECISION_CANDIDATE_MAX_COUNT; ++cand_index)
        context_ptr->fast_candidate_ptr_array[cand_index]->md_rate_estimation_ptr =
            context_ptr->md_rate_estimation_ptr;

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
    context_ptr->qp_index =
        pcs_ptr->parent_pcs_ptr->frm_hdr.delta_q_params.delta_q_present
            ? sb_qp
            : (uint8_t)pcs_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx;

    av1_lambda_assign_md(pcs_ptr, context_ptr);

    return;
}
