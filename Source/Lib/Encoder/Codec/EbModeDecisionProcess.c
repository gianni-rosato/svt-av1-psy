/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>

#include "EbUtility.h"
#include "EbModeDecisionProcess.h"
#include "EbLambdaRateTables.h"

static void mode_decision_context_dctor(EbPtr p) {
    ModeDecisionContext *obj = (ModeDecisionContext *)p;
    for (int cd = 0; cd < MAX_PAL_CAND; cd++)
        if (obj->palette_cand_array[cd].color_idx_map)
            EB_FREE_ARRAY(obj->palette_cand_array[cd].color_idx_map);
    for (uint32_t cand_index = 0; cand_index < MODE_DECISION_CANDIDATE_MAX_COUNT; ++cand_index) {
        if (obj->fast_candidate_ptr_array[cand_index]->palette_info.color_idx_map)
            for (uint32_t coded_leaf_index = 0; coded_leaf_index < BLOCK_MAX_COUNT_SB_128;
                 ++coded_leaf_index)
                if (obj->md_blk_arr_nsq[coded_leaf_index].palette_info.color_idx_map)
                    EB_FREE_ARRAY(obj->md_blk_arr_nsq[coded_leaf_index].palette_info.color_idx_map);
        EB_FREE_ARRAY(obj->fast_candidate_ptr_array[cand_index]->palette_info.color_idx_map);
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
    if (obj->hbd_mode_decision > EB_8_BIT_MD) EB_FREE_ALIGNED_ARRAY(obj->cfl_temp_luma_recon16bit);
    if (obj->hbd_mode_decision != EB_10_BIT_MD) EB_FREE_ALIGNED_ARRAY(obj->cfl_temp_luma_recon);
    if (obj->is_md_rate_estimation_ptr_owner) EB_FREE_ARRAY(obj->md_rate_estimation_ptr);
    EB_FREE_ARRAY(obj->fast_candidate_array);
    EB_FREE_ARRAY(obj->fast_candidate_ptr_array);
    EB_FREE_ARRAY(obj->fast_cost_array);
    EB_FREE_ARRAY(obj->full_cost_array);
    EB_FREE_ARRAY(obj->full_cost_skip_ptr);
    EB_FREE_ARRAY(obj->full_cost_merge_ptr);
    if (obj->md_local_blk_unit) {
        if (obj->hbd_mode_decision > EB_8_BIT_MD) {
            EB_FREE_ARRAY(obj->md_local_blk_unit[0].neigh_left_recon_16bit[0]);
            EB_FREE_ARRAY(obj->md_local_blk_unit[0].neigh_top_recon_16bit[0]);
        }
        if (obj->hbd_mode_decision != EB_10_BIT_MD) {
            EB_FREE_ARRAY(obj->md_local_blk_unit[0].neigh_left_recon[0]);
            EB_FREE_ARRAY(obj->md_local_blk_unit[0].neigh_top_recon[0]);
        }
    }
    if (obj->md_blk_arr_nsq) {
        EB_FREE_ARRAY(obj->md_blk_arr_nsq[0].av1xd);
    }
    EB_FREE_ARRAY(obj->md_local_blk_unit);
    EB_FREE_ARRAY(obj->md_blk_arr_nsq);
    EB_FREE_ARRAY(obj->md_ep_pipe_sb);
#if DEPTH_PART_CLEAN_UP
    EB_FREE_ARRAY(obj->mdc_sb_array);
#endif
#if UNIFY_TXT
    for (uint32_t txt_itr = 0; txt_itr < TX_TYPES; ++txt_itr) {
        EB_DELETE(obj->recon_coeff_ptr[txt_itr]);
        EB_DELETE(obj->recon_ptr[txt_itr]);
    }
#endif
#if CAND_MEM_OPT
    EB_DELETE(obj->prediction_ptr_temp);
    EB_DELETE(obj->cfl_temp_prediction_ptr);
    EB_DELETE(obj->residual_quant_coeff_ptr);
#endif

#if MEM_OPT_MD_BUF_DESC
    EB_DELETE(obj->temp_residual_ptr);
    EB_DELETE(obj->temp_recon_ptr);
#endif
}

/******************************************************
 * Mode Decision Context Constructor
 ******************************************************/
EbErrorType mode_decision_context_ctor(ModeDecisionContext *context_ptr, EbColorFormat color_format,
                                       EbFifo *mode_decision_configuration_input_fifo_ptr,
                                       EbFifo *mode_decision_output_fifo_ptr,
                                       uint8_t enable_hbd_mode_decision, uint8_t cfg_palette) {
    uint32_t buffer_index;
    uint32_t cand_index;

    (void)color_format;

    context_ptr->dctor             = mode_decision_context_dctor;
    context_ptr->hbd_mode_decision = enable_hbd_mode_decision;

    // Input/Output System Resource Manager FIFOs
    context_ptr->mode_decision_configuration_input_fifo_ptr =
        mode_decision_configuration_input_fifo_ptr;
    context_ptr->mode_decision_output_fifo_ptr = mode_decision_output_fifo_ptr;

    // Cfl scratch memory
    if (context_ptr->hbd_mode_decision > EB_8_BIT_MD)
        EB_MALLOC_ALIGNED(context_ptr->cfl_temp_luma_recon16bit, sizeof(uint16_t) * 128 * 128);
    if (context_ptr->hbd_mode_decision != EB_10_BIT_MD)
        EB_MALLOC_ALIGNED(context_ptr->cfl_temp_luma_recon, sizeof(uint8_t) * 128 * 128);
    // MD rate Estimation tables
    EB_MALLOC_ARRAY(context_ptr->md_rate_estimation_ptr, 1);
    context_ptr->is_md_rate_estimation_ptr_owner = EB_TRUE;

    EB_MALLOC_ARRAY(context_ptr->md_local_blk_unit, BLOCK_MAX_COUNT_SB_128);
    EB_MALLOC_ARRAY(context_ptr->md_blk_arr_nsq, BLOCK_MAX_COUNT_SB_128);
    EB_MALLOC_ARRAY(context_ptr->md_ep_pipe_sb, BLOCK_MAX_COUNT_SB_128);

    // Fast Candidate Array
    EB_MALLOC_ARRAY(context_ptr->fast_candidate_array, MODE_DECISION_CANDIDATE_MAX_COUNT);

    EB_MALLOC_ARRAY(context_ptr->fast_candidate_ptr_array, MODE_DECISION_CANDIDATE_MAX_COUNT);

    for (cand_index = 0; cand_index < MODE_DECISION_CANDIDATE_MAX_COUNT; ++cand_index) {
        context_ptr->fast_candidate_ptr_array[cand_index] =
            &context_ptr->fast_candidate_array[cand_index];
        context_ptr->fast_candidate_ptr_array[cand_index]->md_rate_estimation_ptr =
            context_ptr->md_rate_estimation_ptr;
        if (cfg_palette)
            EB_MALLOC_ARRAY(
                context_ptr->fast_candidate_ptr_array[cand_index]->palette_info.color_idx_map,
                MAX_PALETTE_SQUARE);
        else
            context_ptr->fast_candidate_ptr_array[cand_index]->palette_info.color_idx_map = NULL;
    }
    for (int cd = 0; cd < MAX_PAL_CAND; cd++)
        if (cfg_palette)
            EB_MALLOC_ARRAY(context_ptr->palette_cand_array[cd].color_idx_map, MAX_PALETTE_SQUARE);
        else
            context_ptr->palette_cand_array[cd].color_idx_map = NULL;
    // Transform and Quantization Buffers
    EB_NEW(context_ptr->trans_quant_buffers_ptr, eb_trans_quant_buffers_ctor);

    // Cost Arrays
    EB_MALLOC_ARRAY(context_ptr->fast_cost_array, MAX_NFL_BUFF);
    EB_MALLOC_ARRAY(context_ptr->full_cost_array, MAX_NFL_BUFF);
    EB_MALLOC_ARRAY(context_ptr->full_cost_skip_ptr, MAX_NFL_BUFF);
    EB_MALLOC_ARRAY(context_ptr->full_cost_merge_ptr, MAX_NFL_BUFF);
    // Candidate Buffers
    EB_ALLOC_PTR_ARRAY(context_ptr->candidate_buffer_ptr_array, MAX_NFL_BUFF);
    for (buffer_index = 0; buffer_index < MAX_NFL_BUFF; ++buffer_index) {
        EB_NEW(context_ptr->candidate_buffer_ptr_array[buffer_index],
               mode_decision_candidate_buffer_ctor,
               context_ptr->hbd_mode_decision ? EB_10BIT : EB_8BIT,
               &(context_ptr->fast_cost_array[buffer_index]),
               &(context_ptr->full_cost_array[buffer_index]),
               &(context_ptr->full_cost_skip_ptr[buffer_index]),
               &(context_ptr->full_cost_merge_ptr[buffer_index]));
    }
    EB_NEW(context_ptr->candidate_buffer_tx_depth_1,
           mode_decision_scratch_candidate_buffer_ctor,
           context_ptr->hbd_mode_decision ? EB_10BIT : EB_8BIT);

    EB_ALLOC_PTR_ARRAY(context_ptr->candidate_buffer_tx_depth_1->candidate_ptr, 1);

    EB_NEW(context_ptr->candidate_buffer_tx_depth_2,
           mode_decision_scratch_candidate_buffer_ctor,
           context_ptr->hbd_mode_decision ? EB_10BIT : EB_8BIT);

    EB_ALLOC_PTR_ARRAY(context_ptr->candidate_buffer_tx_depth_2->candidate_ptr, 1);
    context_ptr->md_local_blk_unit[0].neigh_left_recon[0]       = NULL;
    context_ptr->md_local_blk_unit[0].neigh_top_recon[0]        = NULL;
    context_ptr->md_local_blk_unit[0].neigh_left_recon_16bit[0] = NULL;
    context_ptr->md_local_blk_unit[0].neigh_top_recon_16bit[0]  = NULL;
    uint16_t sz                                                 = sizeof(uint16_t);
    if (context_ptr->hbd_mode_decision > EB_8_BIT_MD) {
        EB_MALLOC_ARRAY(context_ptr->md_local_blk_unit[0].neigh_left_recon_16bit[0],
                        BLOCK_MAX_COUNT_SB_128 * 128 * 3 * sz);
        EB_MALLOC_ARRAY(context_ptr->md_local_blk_unit[0].neigh_top_recon_16bit[0],
                        BLOCK_MAX_COUNT_SB_128 * 128 * 3 * sz);
    }
    if (context_ptr->hbd_mode_decision != EB_10_BIT_MD) {
        EB_MALLOC_ARRAY(context_ptr->md_local_blk_unit[0].neigh_left_recon[0],
                        BLOCK_MAX_COUNT_SB_128 * 128 * 3);
        EB_MALLOC_ARRAY(context_ptr->md_local_blk_unit[0].neigh_top_recon[0],
                        BLOCK_MAX_COUNT_SB_128 * 128 * 3);
    }
    uint32_t coded_leaf_index;
    for (coded_leaf_index = 0; coded_leaf_index < BLOCK_MAX_COUNT_SB_128; ++coded_leaf_index) {
        for (int i = 0; i < 3; i++) {
            size_t offset = (coded_leaf_index * 128 * 3 + i * 128) * sz;
            context_ptr->md_local_blk_unit[coded_leaf_index].neigh_left_recon_16bit[i] =
                context_ptr->md_local_blk_unit[0].neigh_left_recon_16bit[0] + offset;
            context_ptr->md_local_blk_unit[coded_leaf_index].neigh_top_recon_16bit[i] =
                context_ptr->md_local_blk_unit[0].neigh_top_recon_16bit[0] + offset;
        }
        for (int i = 0; i < 3; i++) {
            size_t offset = coded_leaf_index * 128 * 3 + i * 128;
            context_ptr->md_local_blk_unit[coded_leaf_index].neigh_left_recon[i] =
                context_ptr->md_local_blk_unit[0].neigh_left_recon[0] + offset;
            context_ptr->md_local_blk_unit[coded_leaf_index].neigh_top_recon[i] =
                context_ptr->md_local_blk_unit[0].neigh_top_recon[0] + offset;
        }
    }
    context_ptr->md_blk_arr_nsq[0].av1xd                     = NULL;
    EB_MALLOC_ARRAY(context_ptr->md_blk_arr_nsq[0].av1xd, BLOCK_MAX_COUNT_SB_128);

#if DEPTH_PART_CLEAN_UP
    EB_MALLOC_ARRAY(context_ptr->mdc_sb_array, 1);
#endif
    for (coded_leaf_index = 0; coded_leaf_index < BLOCK_MAX_COUNT_SB_128; ++coded_leaf_index) {
        context_ptr->md_blk_arr_nsq[coded_leaf_index].av1xd =
            context_ptr->md_blk_arr_nsq[0].av1xd + coded_leaf_index;
        context_ptr->md_blk_arr_nsq[coded_leaf_index].segment_id = 0;
        if (cfg_palette)
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
                   eb_picture_buffer_desc_ctor,
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
                   eb_picture_buffer_desc_ctor,
                   (EbPtr)&init_data);
        }
#endif
    }
    EB_MALLOC_ARRAY(context_ptr->ref_best_cost_sq_table, MAX_REF_TYPE_CAND);
    EB_MALLOC_ARRAY(context_ptr->ref_best_ref_sq_table, MAX_REF_TYPE_CAND);
    EB_MALLOC_ARRAY(context_ptr->above_txfm_context, (MAX_SB_SIZE >> MI_SIZE_LOG2));
    EB_MALLOC_ARRAY(context_ptr->left_txfm_context, (MAX_SB_SIZE >> MI_SIZE_LOG2));
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

extern void lambda_assign_low_delay(uint32_t *fast_lambda, uint32_t *full_lambda,
                                    uint32_t *fast_chroma_lambda, uint32_t *full_chroma_lambda,
                                    uint32_t *full_chroma_lambda_sao,
                                    uint8_t qp_hierarchical_position, uint8_t qp, uint8_t chroma_qp)

{
    if (qp_hierarchical_position == 0) {
        *fast_lambda            = lambda_mode_decision_ld_sad[qp];
        *fast_chroma_lambda     = lambda_mode_decision_ld_sad[qp];
        *full_lambda            = lambda_mode_decision_ld_sse[qp];
        *full_chroma_lambda     = lambda_mode_decision_ld_sse[qp];
        *full_chroma_lambda_sao = lambda_mode_decision_ld_sse[chroma_qp];
    } else { // Hierarchical postions 1, 2, 3, 4, 5
        *fast_lambda            = lambda_mode_decision_ld_sad_qp_scaling[qp];
        *fast_chroma_lambda     = lambda_mode_decision_ld_sad_qp_scaling[qp];
        *full_lambda            = lambda_mode_decision_ld_sse_qp_scaling[qp];
        *full_chroma_lambda     = lambda_mode_decision_ld_sse_qp_scaling[qp];
        *full_chroma_lambda_sao = lambda_mode_decision_ld_sse_qp_scaling[chroma_qp];
    }
}

void lambda_assign_random_access(uint32_t *fast_lambda, uint32_t *full_lambda,
                                 uint32_t *fast_chroma_lambda, uint32_t *full_chroma_lambda,
                                 uint32_t *full_chroma_lambda_sao, uint8_t qp_hierarchical_position,
                                 uint8_t qp, uint8_t chroma_qp)

{
    if (qp_hierarchical_position == 0) {
        *fast_lambda            = lambda_mode_decision_ra_sad[qp];
        *fast_chroma_lambda     = lambda_mode_decision_ra_sad[qp];
        *full_lambda            = lambda_mode_decision_ra_sse[qp];
        *full_chroma_lambda     = lambda_mode_decision_ra_sse[qp];
        *full_chroma_lambda_sao = lambda_mode_decision_ra_sse[chroma_qp];
    } else if (qp_hierarchical_position < 3) { // Hierarchical postions 1, 2

        *fast_lambda            = lambda_mode_decision_ra_sad_qp_scaling_l1[qp];
        *fast_chroma_lambda     = lambda_mode_decision_ra_sad_qp_scaling_l1[qp];
        *full_lambda            = lambda_mode_decision_ra_sse_qp_scaling_l1[qp];
        *full_chroma_lambda     = lambda_mode_decision_ra_sse_qp_scaling_l1[qp];
        *full_chroma_lambda_sao = lambda_mode_decision_ra_sse_qp_scaling_l1[chroma_qp];
    } else { // Hierarchical postions 3, 4, 5
        *fast_lambda            = lambda_mode_decision_ra_sad_qp_scaling_l3[qp];
        *fast_chroma_lambda     = lambda_mode_decision_ra_sad_qp_scaling_l3[qp];
        *full_lambda            = lambda_mode_decision_ra_sse_qp_scaling_l3[qp];
        *full_chroma_lambda     = lambda_mode_decision_ra_sse_qp_scaling_l3[qp];
        *full_chroma_lambda_sao = lambda_mode_decision_ra_sse_qp_scaling_l3[chroma_qp];
    }
}

void lambda_assign_i_slice(uint32_t *fast_lambda, uint32_t *full_lambda,
                           uint32_t *fast_chroma_lambda, uint32_t *full_chroma_lambda,
                           uint32_t *full_chroma_lambda_sao, uint8_t qp_hierarchical_position,
                           uint8_t qp, uint8_t chroma_qp)

{
    if (qp_hierarchical_position == 0) {
        *fast_lambda            = lambda_mode_decision_i_slice_sad[qp];
        *fast_chroma_lambda     = lambda_mode_decision_i_slice_sad[qp];
        *full_lambda            = lambda_mode_decision_i_slice_sse[qp];
        *full_chroma_lambda     = lambda_mode_decision_i_slice_sse[qp];
        *full_chroma_lambda_sao = lambda_mode_decision_i_slice_sse[chroma_qp];
    }
}
const EbLambdaAssignFunc lambda_assignment_function_table[4] = {
    lambda_assign_low_delay, // low delay P
    lambda_assign_low_delay, // low delay b
    lambda_assign_random_access, // Random Access
    lambda_assign_i_slice // I_SLICE
};

void av1_lambda_assign_md(
    ModeDecisionContext   *context_ptr)
{
        context_ptr->full_lambda_md[0] = av1_lambda_mode_decision8_bit_sse[context_ptr->qp_index];
        context_ptr->fast_lambda_md[0] = av1_lambda_mode_decision8_bit_sad[context_ptr->qp_index];

        context_ptr->full_lambda_md[1] = av1lambda_mode_decision10_bit_sse[context_ptr->qp_index];
        context_ptr->fast_lambda_md[1] = av1lambda_mode_decision10_bit_sad[context_ptr->qp_index];

        context_ptr->full_lambda_md[1] *= 16;
        context_ptr->fast_lambda_md[1] *= 4;
}
void av1_lambda_assign(uint32_t *fast_lambda, uint32_t *full_lambda, uint8_t bit_depth, uint16_t qp_index,
                       EbBool multiply_lambda) {
    if (bit_depth == 8) {
        *full_lambda = av1_lambda_mode_decision8_bit_sse[qp_index];
        *fast_lambda = av1_lambda_mode_decision8_bit_sad[qp_index];
    } else if (bit_depth == 10) {
        *full_lambda = av1lambda_mode_decision10_bit_sse[qp_index];
        *fast_lambda = av1lambda_mode_decision10_bit_sad[qp_index];
        if (multiply_lambda) {
            *full_lambda *= 16;
            *fast_lambda *= 4;
        }
    } else if (bit_depth == 12) {
        *full_lambda = av1lambda_mode_decision12_bit_sse[qp_index];
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
    uint16_t picture_qp   = pcs_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx;
    context_ptr->qp       = picture_qp;
    context_ptr->qp_index = context_ptr->qp;
    // Asuming cb and cr offset to be the same for chroma QP in both slice and pps for lambda computation
    context_ptr->chroma_qp = (uint8_t)context_ptr->qp;
    context_ptr->qp_index  = (uint8_t)frm_hdr->quantization_params.base_q_idx;
    av1_lambda_assign_md(context_ptr);
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
    context_ptr->coeff_est_entropy_coder_ptr = pcs_ptr->coeff_est_entropy_coder_ptr;

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
    (void)pcs_ptr;
    //Disable Lambda update per SB
    context_ptr->qp = sb_qp;
    // Asuming cb and cr offset to be the same for chroma QP in both slice and pps for lambda computation

    context_ptr->chroma_qp = (uint8_t)context_ptr->qp;

    /* Note(CHKN) : when Qp modulation varies QP on a sub-SB(CU) basis,  Lamda has to change based on Cu->QP , and then this code has to move inside the CU loop in MD */

    // Lambda Assignement
    context_ptr->qp_index =
        pcs_ptr->parent_pcs_ptr->frm_hdr.delta_q_params.delta_q_present
            ? (uint8_t)quantizer_to_qindex[sb_qp]
            : (uint8_t)pcs_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx;

    av1_lambda_assign_md(context_ptr);

    return;
}
