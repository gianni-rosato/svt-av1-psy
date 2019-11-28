/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>

#include "EbUtility.h"
#include "EbModeDecisionProcess.h"
#include "EbLambdaRateTables.h"

static void mode_decision_context_dctor(EbPtr p)
{
    ModeDecisionContext* obj = (ModeDecisionContext*)p;
#if PAL_SUP
    for (int cd = 0; cd < MAX_PAL_CAND; cd++)
        if (obj->palette_cand_array[cd].color_idx_map)
            EB_FREE_ARRAY(obj->palette_cand_array[cd].color_idx_map);
    for (uint32_t candidateIndex = 0; candidateIndex < MODE_DECISION_CANDIDATE_MAX_COUNT; ++candidateIndex)
        if (obj->fast_candidate_ptr_array[candidateIndex]->palette_info.color_idx_map)
            EB_FREE_ARRAY(obj->fast_candidate_ptr_array[candidateIndex]->palette_info.color_idx_map);
    for (uint32_t codedLeafIndex = 0; codedLeafIndex < BLOCK_MAX_COUNT_SB_128; ++codedLeafIndex)
        if (obj->md_cu_arr_nsq[codedLeafIndex].palette_info.color_idx_map)
            EB_FREE_ARRAY(obj->md_cu_arr_nsq[codedLeafIndex].palette_info.color_idx_map);
#endif
    EB_FREE_ARRAY(obj->ref_best_ref_sq_table);
    EB_FREE_ARRAY(obj->ref_best_cost_sq_table);
#if ENHANCE_ATB
    EB_FREE_ARRAY(obj->above_txfm_context);
    EB_FREE_ARRAY(obj->left_txfm_context);
#endif
#if NO_ENCDEC //SB128_TODO to upgrade
    int codedLeafIndex;
    for (codedLeafIndex = 0; codedLeafIndex < BLOCK_MAX_COUNT_SB_128; ++codedLeafIndex) {
        EB_DELETE(obj->md_cu_arr_nsq[codedLeafIndex].recon_tmp);
        EB_DELETE(obj->md_cu_arr_nsq[codedLeafIndex].coeff_tmp);

   }
#endif
    EB_DELETE_PTR_ARRAY(obj->candidate_buffer_ptr_array, MAX_NFL_BUFF);

#if ENHANCE_ATB
    EB_FREE_ARRAY(obj->scratch_candidate_buffer->candidate_ptr);
    EB_DELETE(obj->scratch_candidate_buffer);
#endif

    EB_DELETE(obj->trans_quant_buffers_ptr);
#if !HBD_CLEAN_UP // cfl_temp_luma_recon16bit
    if (obj->hbd_mode_decision)
#else
    if (obj->hbd_mode_decision > EB_8_BIT_MD)
#endif
        EB_FREE_ALIGNED_ARRAY(obj->cfl_temp_luma_recon16bit);
#if !HBD_CLEAN_UP
    else
#else
    if (obj->hbd_mode_decision != EB_10_BIT_MD)
#endif
        EB_FREE_ALIGNED_ARRAY(obj->cfl_temp_luma_recon);
    EB_FREE(obj->transform_inner_array_ptr);
    if (obj->is_md_rate_estimation_ptr_owner)
        EB_FREE_ARRAY(obj->md_rate_estimation_ptr);
    EB_FREE_ARRAY(obj->fast_candidate_array);
    EB_FREE_ARRAY(obj->fast_candidate_ptr_array);
    EB_FREE_ARRAY(obj->fast_cost_array);
    EB_FREE_ARRAY(obj->full_cost_array);
    EB_FREE_ARRAY(obj->full_cost_skip_ptr);
    EB_FREE_ARRAY(obj->full_cost_merge_ptr);
    if (obj->md_cu_arr_nsq) {
        EB_FREE_ARRAY(obj->md_cu_arr_nsq[0].av1xd);
#if !HBD_CLEAN_UP // neigh_left_recon_16bit neigh_top_recon_16bit
        if (obj->hbd_mode_decision) {
#else
        if (obj->hbd_mode_decision > EB_8_BIT_MD) {
#endif
            EB_FREE_ARRAY(obj->md_cu_arr_nsq[0].neigh_left_recon_16bit[0]);
            EB_FREE_ARRAY(obj->md_cu_arr_nsq[0].neigh_top_recon_16bit[0]);
        }
#if HBD_CLEAN_UP
        if (obj->hbd_mode_decision != EB_10_BIT_MD) {
#else
        else {

#endif
            EB_FREE_ARRAY(obj->md_cu_arr_nsq[0].neigh_left_recon[0]);
            EB_FREE_ARRAY(obj->md_cu_arr_nsq[0].neigh_top_recon[0]);
        }
    }

    EB_FREE_ARRAY(obj->md_local_cu_unit);
    EB_FREE_ARRAY(obj->md_cu_arr_nsq);
    EB_FREE_ARRAY(obj->md_ep_pipe_sb);
}

/******************************************************
 * Mode Decision Context Constructor
 ******************************************************/
EbErrorType mode_decision_context_ctor(
    ModeDecisionContext  *context_ptr,
    EbColorFormat         color_format,
    EbFifo                *mode_decision_configuration_input_fifo_ptr,
    EbFifo                *mode_decision_output_fifo_ptr,
    uint8_t                enable_hbd_mode_decision
#if PAL_SUP
    ,uint8_t                 cfg_palette
#endif
)
{
    uint32_t bufferIndex;
    uint32_t candidateIndex;

    (void)color_format;

    context_ptr->dctor = mode_decision_context_dctor;
    context_ptr->hbd_mode_decision = enable_hbd_mode_decision;

    // Input/Output System Resource Manager FIFOs
    context_ptr->mode_decision_configuration_input_fifo_ptr = mode_decision_configuration_input_fifo_ptr;
    context_ptr->mode_decision_output_fifo_ptr = mode_decision_output_fifo_ptr;

    // Trasform Scratch Memory
    EB_MALLOC(context_ptr->transform_inner_array_ptr, 3120); //refer to EbInvTransform_SSE2.as. case 32x32

    // Cfl scratch memory
#if !HBD_CLEAN_UP // cfl_temp_luma_recon16bit
    if (context_ptr->hbd_mode_decision) {
#else
    if (context_ptr->hbd_mode_decision > EB_8_BIT_MD)
#endif
        EB_MALLOC_ALIGNED(context_ptr->cfl_temp_luma_recon16bit, sizeof(uint16_t) * 128 * 128);
#if !HBD_CLEAN_UP
    } else {
#else
    if (context_ptr->hbd_mode_decision != EB_10_BIT_MD)
#endif
        EB_MALLOC_ALIGNED(context_ptr->cfl_temp_luma_recon, sizeof(uint8_t) * 128 * 128);
#if !HBD_CLEAN_UP
    }
#endif
    // MD rate Estimation tables
    EB_MALLOC_ARRAY(context_ptr->md_rate_estimation_ptr, 1);
    context_ptr->is_md_rate_estimation_ptr_owner = EB_TRUE;

    EB_MALLOC_ARRAY(context_ptr->md_local_cu_unit, BLOCK_MAX_COUNT_SB_128);
    EB_MALLOC_ARRAY(context_ptr->md_cu_arr_nsq, BLOCK_MAX_COUNT_SB_128);
    EB_MALLOC_ARRAY(context_ptr->md_ep_pipe_sb, BLOCK_MAX_COUNT_SB_128);

    // Fast Candidate Array
    EB_MALLOC_ARRAY(context_ptr->fast_candidate_array, MODE_DECISION_CANDIDATE_MAX_COUNT);

    EB_MALLOC_ARRAY(context_ptr->fast_candidate_ptr_array, MODE_DECISION_CANDIDATE_MAX_COUNT);

    for (candidateIndex = 0; candidateIndex < MODE_DECISION_CANDIDATE_MAX_COUNT; ++candidateIndex) {
        context_ptr->fast_candidate_ptr_array[candidateIndex] = &context_ptr->fast_candidate_array[candidateIndex];
        context_ptr->fast_candidate_ptr_array[candidateIndex]->md_rate_estimation_ptr = context_ptr->md_rate_estimation_ptr;
#if PAL_SUP
        if (cfg_palette)
            EB_MALLOC_ARRAY(context_ptr->fast_candidate_ptr_array[candidateIndex]->palette_info.color_idx_map, MAX_PALETTE_SQUARE);
        else
            context_ptr->fast_candidate_ptr_array[candidateIndex]->palette_info.color_idx_map = NULL;
#endif
    }

#if PAL_SUP
    for (int cd = 0; cd < MAX_PAL_CAND; cd++)
        if (cfg_palette)
            EB_MALLOC_ARRAY(context_ptr->palette_cand_array[cd].color_idx_map, MAX_PALETTE_SQUARE);
        else
            context_ptr->palette_cand_array[cd].color_idx_map = NULL;
#endif
    // Transform and Quantization Buffers
    EB_NEW(
        context_ptr->trans_quant_buffers_ptr,
        eb_trans_quant_buffers_ctor);

    // Cost Arrays
    EB_MALLOC_ARRAY(context_ptr->fast_cost_array, MAX_NFL_BUFF);
    EB_MALLOC_ARRAY(context_ptr->full_cost_array, MAX_NFL_BUFF);
    EB_MALLOC_ARRAY(context_ptr->full_cost_skip_ptr, MAX_NFL_BUFF);
    EB_MALLOC_ARRAY(context_ptr->full_cost_merge_ptr, MAX_NFL_BUFF);
    // Candidate Buffers
    EB_ALLOC_PTR_ARRAY(context_ptr->candidate_buffer_ptr_array, MAX_NFL_BUFF);
    for (bufferIndex = 0; bufferIndex < MAX_NFL_BUFF; ++bufferIndex) {

        EB_NEW(
            context_ptr->candidate_buffer_ptr_array[bufferIndex],
            mode_decision_candidate_buffer_ctor,
            context_ptr->hbd_mode_decision ? EB_10BIT : EB_8BIT,
            &(context_ptr->fast_cost_array[bufferIndex]),
            &(context_ptr->full_cost_array[bufferIndex]),
            &(context_ptr->full_cost_skip_ptr[bufferIndex]),
            &(context_ptr->full_cost_merge_ptr[bufferIndex])
        );
    }
#if ENHANCE_ATB
    EB_NEW(
        context_ptr->scratch_candidate_buffer,
        mode_decision_scratch_candidate_buffer_ctor,
        context_ptr->hbd_mode_decision ? EB_10BIT : EB_8BIT);

    EB_ALLOC_PTR_ARRAY(context_ptr->scratch_candidate_buffer->candidate_ptr, 1);
#endif
    context_ptr->md_cu_arr_nsq[0].av1xd = NULL;
    context_ptr->md_cu_arr_nsq[0].neigh_left_recon[0] = NULL;
    context_ptr->md_cu_arr_nsq[0].neigh_top_recon[0] = NULL;
    context_ptr->md_cu_arr_nsq[0].neigh_left_recon_16bit[0] = NULL;
    context_ptr->md_cu_arr_nsq[0].neigh_top_recon_16bit[0] = NULL;
    EB_MALLOC_ARRAY(context_ptr->md_cu_arr_nsq[0].av1xd, BLOCK_MAX_COUNT_SB_128);
    uint16_t sz = sizeof(uint16_t);
#if !HBD_CLEAN_UP // neigh_left_recon_16bit neigh_top_recon_16bit
    if (context_ptr->hbd_mode_decision) {
#else
    if (context_ptr->hbd_mode_decision > EB_8_BIT_MD){
#endif
        EB_MALLOC_ARRAY(context_ptr->md_cu_arr_nsq[0].neigh_left_recon_16bit[0], BLOCK_MAX_COUNT_SB_128 * 128 * 3 * sz);
        EB_MALLOC_ARRAY(context_ptr->md_cu_arr_nsq[0].neigh_top_recon_16bit[0], BLOCK_MAX_COUNT_SB_128 * 128 * 3 * sz);
    }
#if HBD_CLEAN_UP
    if (context_ptr->hbd_mode_decision != EB_10_BIT_MD){
#else
    else {
#endif
        EB_MALLOC_ARRAY(context_ptr->md_cu_arr_nsq[0].neigh_left_recon[0], BLOCK_MAX_COUNT_SB_128 * 128 * 3);
        EB_MALLOC_ARRAY(context_ptr->md_cu_arr_nsq[0].neigh_top_recon[0], BLOCK_MAX_COUNT_SB_128 * 128 * 3);
    }
    uint32_t codedLeafIndex, tu_index;
    for (codedLeafIndex = 0; codedLeafIndex < BLOCK_MAX_COUNT_SB_128; ++codedLeafIndex) {
        for (tu_index = 0; tu_index < TRANSFORM_UNIT_MAX_COUNT; ++tu_index)
            context_ptr->md_cu_arr_nsq[codedLeafIndex].transform_unit_array[tu_index].tu_index = tu_index;
        const BlockGeom * blk_geom = get_blk_geom_mds(codedLeafIndex);
        UNUSED(blk_geom);
        context_ptr->md_cu_arr_nsq[codedLeafIndex].av1xd = context_ptr->md_cu_arr_nsq[0].av1xd + codedLeafIndex;
#if !HBD_CLEAN_UP
        if (context_ptr->hbd_mode_decision) {
#endif
             for (int i = 0; i < 3; i++) {
                size_t offset = (codedLeafIndex * 128 * 3 + i * 128) * sz;
                context_ptr->md_cu_arr_nsq[codedLeafIndex].neigh_left_recon_16bit[i] = context_ptr->md_cu_arr_nsq[0].neigh_left_recon_16bit[0] + offset;
                context_ptr->md_cu_arr_nsq[codedLeafIndex].neigh_top_recon_16bit[i] = context_ptr->md_cu_arr_nsq[0].neigh_top_recon_16bit[0] + offset;
            }
#if !HBD_CLEAN_UP
        } else {
#endif
             for (int i = 0; i < 3; i++) {
                size_t offset = codedLeafIndex * 128 * 3 + i * 128;
                context_ptr->md_cu_arr_nsq[codedLeafIndex].neigh_left_recon[i] = context_ptr->md_cu_arr_nsq[0].neigh_left_recon[0] + offset;
                context_ptr->md_cu_arr_nsq[codedLeafIndex].neigh_top_recon[i] = context_ptr->md_cu_arr_nsq[0].neigh_top_recon[0] + offset;
#if !HBD_CLEAN_UP
            }
#endif
        }
#if PAL_SUP
        if (cfg_palette)
            EB_MALLOC_ARRAY(context_ptr->md_cu_arr_nsq[codedLeafIndex].palette_info.color_idx_map, MAX_PALETTE_SQUARE);
        else
            context_ptr->md_cu_arr_nsq[codedLeafIndex].palette_info.color_idx_map = NULL;
#endif
#if NO_ENCDEC //SB128_TODO to upgrade
        {
            EbPictureBufferDescInitData initData;

            initData.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
            initData.max_width = SB_STRIDE_Y;
            initData.max_height = SB_STRIDE_Y;
            initData.bit_depth = EB_32BIT;
            initData.color_format = EB_YUV420;
            initData.left_padding = 0;
            initData.right_padding = 0;
            initData.top_padding = 0;
            initData.bot_padding = 0;
            initData.split_mode = EB_FALSE;

            EB_NEW(
                context_ptr->md_cu_arr_nsq[codedLeafIndex].coeff_tmp,
                eb_picture_buffer_desc_ctor,
                (EbPtr)&initData);

            initData.buffer_enable_mask = PICTURE_BUFFER_DESC_FULL_MASK;
            initData.max_width = SB_STRIDE_Y;
            initData.max_height = SB_STRIDE_Y;
            initData.bit_depth = EB_8BIT;
            initData.color_format = EB_YUV420;
            initData.left_padding = 0;
            initData.right_padding = 0;
            initData.top_padding = 0;
            initData.bot_padding = 0;
            initData.split_mode = EB_FALSE;

            EB_NEW(
                context_ptr->md_cu_arr_nsq[codedLeafIndex].recon_tmp,
                eb_picture_buffer_desc_ctor,
                (EbPtr)&initData);
        }
#endif
    }
    EB_MALLOC_ARRAY(context_ptr->ref_best_cost_sq_table, MAX_REF_TYPE_CAND);
    EB_MALLOC_ARRAY(context_ptr->ref_best_ref_sq_table, MAX_REF_TYPE_CAND);
#if ENHANCE_ATB
    EB_MALLOC_ARRAY(context_ptr->above_txfm_context, (MAX_SB_SIZE >> MI_SIZE_LOG2));
    EB_MALLOC_ARRAY(context_ptr->left_txfm_context, (MAX_SB_SIZE >> MI_SIZE_LOG2));
#endif
    return EB_ErrorNone;
}

/**************************************************
 * Reset Mode Decision Neighbor Arrays
 *************************************************/
void reset_mode_decision_neighbor_arrays(PictureControlSet *picture_control_set_ptr)
{
    uint8_t depth;
    for (depth = 0; depth < NEIGHBOR_ARRAY_TOTAL_COUNT; depth++) {
        neighbor_array_unit_reset(picture_control_set_ptr->md_intra_luma_mode_neighbor_array[depth]);
        neighbor_array_unit_reset(picture_control_set_ptr->md_intra_chroma_mode_neighbor_array[depth]);
        neighbor_array_unit_reset(picture_control_set_ptr->md_mv_neighbor_array[depth]);
        neighbor_array_unit_reset(picture_control_set_ptr->md_skip_flag_neighbor_array[depth]);
        neighbor_array_unit_reset(picture_control_set_ptr->md_mode_type_neighbor_array[depth]);
        neighbor_array_unit_reset(picture_control_set_ptr->md_leaf_depth_neighbor_array[depth]);
        neighbor_array_unit_reset(picture_control_set_ptr->mdleaf_partition_neighbor_array[depth]);

#if !HBD_CLEAN_UP // md_luma_recon_neighbor_array md_tx_depth_1_luma_recon_neighbor_array
        if (!picture_control_set_ptr->hbd_mode_decision) {
#else
        if (picture_control_set_ptr->hbd_mode_decision != EB_10_BIT_MD) {
#endif
            neighbor_array_unit_reset(picture_control_set_ptr->md_luma_recon_neighbor_array[depth]);
            neighbor_array_unit_reset(picture_control_set_ptr->md_tx_depth_1_luma_recon_neighbor_array[depth]);
            neighbor_array_unit_reset(picture_control_set_ptr->md_cb_recon_neighbor_array[depth]);
            neighbor_array_unit_reset(picture_control_set_ptr->md_cr_recon_neighbor_array[depth]);
        }
#if HBD_CLEAN_UP
        if (picture_control_set_ptr->hbd_mode_decision > EB_8_BIT_MD) {
#else
         else {
#endif

            neighbor_array_unit_reset(picture_control_set_ptr->md_luma_recon_neighbor_array16bit[depth]);
            neighbor_array_unit_reset(picture_control_set_ptr->md_tx_depth_1_luma_recon_neighbor_array16bit[depth]);
            neighbor_array_unit_reset(picture_control_set_ptr->md_cb_recon_neighbor_array16bit[depth]);
            neighbor_array_unit_reset(picture_control_set_ptr->md_cr_recon_neighbor_array16bit[depth]);
        }

        neighbor_array_unit_reset(picture_control_set_ptr->md_skip_coeff_neighbor_array[depth]);
        neighbor_array_unit_reset(picture_control_set_ptr->md_luma_dc_sign_level_coeff_neighbor_array[depth]);
        neighbor_array_unit_reset(picture_control_set_ptr->md_tx_depth_1_luma_dc_sign_level_coeff_neighbor_array[depth]);
        neighbor_array_unit_reset(picture_control_set_ptr->md_cb_dc_sign_level_coeff_neighbor_array[depth]);
        neighbor_array_unit_reset(picture_control_set_ptr->md_cr_dc_sign_level_coeff_neighbor_array[depth]);
        neighbor_array_unit_reset(picture_control_set_ptr->md_txfm_context_array[depth]);
        neighbor_array_unit_reset(picture_control_set_ptr->md_inter_pred_dir_neighbor_array[depth]);
        neighbor_array_unit_reset(picture_control_set_ptr->md_ref_frame_type_neighbor_array[depth]);

        neighbor_array_unit_reset32(picture_control_set_ptr->md_interpolation_type_neighbor_array[depth]);
    }

    return;
}

extern void lambda_assign_low_delay(
    uint32_t                    *fast_lambda,
    uint32_t                    *full_lambda,
    uint32_t                    *fast_chroma_lambda,
    uint32_t                    *full_chroma_lambda,
    uint32_t                    *full_chroma_lambda_sao,
    uint8_t                      qp_hierarchical_position,
    uint8_t                      qp,
    uint8_t                      chroma_qp)

{
    if (qp_hierarchical_position == 0) {
        *fast_lambda = lambda_mode_decision_ld_sad[qp];
        *fast_chroma_lambda = lambda_mode_decision_ld_sad[qp];
        *full_lambda = lambda_mode_decision_ld_sse[qp];
        *full_chroma_lambda = lambda_mode_decision_ld_sse[qp];
        *full_chroma_lambda_sao = lambda_mode_decision_ld_sse[chroma_qp];
    }
    else { // Hierarchical postions 1, 2, 3, 4, 5
        *fast_lambda = lambda_mode_decision_ld_sad_qp_scaling[qp];
        *fast_chroma_lambda = lambda_mode_decision_ld_sad_qp_scaling[qp];
        *full_lambda = lambda_mode_decision_ld_sse_qp_scaling[qp];
        *full_chroma_lambda = lambda_mode_decision_ld_sse_qp_scaling[qp];
        *full_chroma_lambda_sao = lambda_mode_decision_ld_sse_qp_scaling[chroma_qp];
    }
}

void lambda_assign_random_access(
    uint32_t                    *fast_lambda,
    uint32_t                    *full_lambda,
    uint32_t                    *fast_chroma_lambda,
    uint32_t                    *full_chroma_lambda,
    uint32_t                    *full_chroma_lambda_sao,
    uint8_t                      qp_hierarchical_position,
    uint8_t                      qp,
    uint8_t                      chroma_qp)

{
    if (qp_hierarchical_position == 0) {
        *fast_lambda = lambda_mode_decision_ra_sad[qp];
        *fast_chroma_lambda = lambda_mode_decision_ra_sad[qp];
        *full_lambda = lambda_mode_decision_ra_sse[qp];
        *full_chroma_lambda = lambda_mode_decision_ra_sse[qp];
        *full_chroma_lambda_sao = lambda_mode_decision_ra_sse[chroma_qp];
    }
    else if (qp_hierarchical_position < 3) { // Hierarchical postions 1, 2

        *fast_lambda = lambda_mode_decision_ra_sad_qp_scaling_l1[qp];
        *fast_chroma_lambda = lambda_mode_decision_ra_sad_qp_scaling_l1[qp];
        *full_lambda = lambda_mode_decision_ra_sse_qp_scaling_l1[qp];
        *full_chroma_lambda = lambda_mode_decision_ra_sse_qp_scaling_l1[qp];
        *full_chroma_lambda_sao = lambda_mode_decision_ra_sse_qp_scaling_l1[chroma_qp];
    }
    else { // Hierarchical postions 3, 4, 5
        *fast_lambda = lambda_mode_decision_ra_sad_qp_scaling_l3[qp];
        *fast_chroma_lambda = lambda_mode_decision_ra_sad_qp_scaling_l3[qp];
        *full_lambda = lambda_mode_decision_ra_sse_qp_scaling_l3[qp];
        *full_chroma_lambda = lambda_mode_decision_ra_sse_qp_scaling_l3[qp];
        *full_chroma_lambda_sao = lambda_mode_decision_ra_sse_qp_scaling_l3[chroma_qp];
    }
}

void lambdaAssignISlice(
    uint32_t                    *fast_lambda,
    uint32_t                    *full_lambda,
    uint32_t                    *fast_chroma_lambda,
    uint32_t                    *full_chroma_lambda,
    uint32_t                    *full_chroma_lambda_sao,
    uint8_t                      qp_hierarchical_position,
    uint8_t                      qp,
    uint8_t                      chroma_qp)

{
    if (qp_hierarchical_position == 0) {
        *fast_lambda = lambda_mode_decision_i_slice_sad[qp];
        *fast_chroma_lambda = lambda_mode_decision_i_slice_sad[qp];
        *full_lambda = lambda_mode_decision_i_slice_sse[qp];
        *full_chroma_lambda = lambda_mode_decision_i_slice_sse[qp];
        *full_chroma_lambda_sao = lambda_mode_decision_i_slice_sse[chroma_qp];
    }
}
const EbLambdaAssignFunc lambda_assignment_function_table[4] = {
    lambda_assign_low_delay, // low delay P
    lambda_assign_low_delay, // low delay B
    lambda_assign_random_access, // Random Access
    lambdaAssignISlice // I_SLICE
};

void Av1lambdaAssign(
    uint32_t                    *fast_lambda,
    uint32_t                    *full_lambda,
    uint32_t                    *fast_chroma_lambda,
    uint32_t                    *full_chroma_lambda,
    uint8_t                      bit_depth,
    uint16_t                     qp_index,
    EbBool                       hbd_mode_decision)
{
    if (bit_depth == 8) {
        *full_lambda = av1_lambda_mode_decision8_bit_sse[qp_index];
        *fast_lambda = av1_lambda_mode_decision8_bit_sad[qp_index];
    }
    else if (bit_depth == 10) {
        *full_lambda = av1lambda_mode_decision10_bit_sse[qp_index];
        *fast_lambda = av1lambda_mode_decision10_bit_sad[qp_index];
        if (hbd_mode_decision) {
            *full_lambda *= 16;
            *fast_lambda *= 4;
        }
    }
    else if (bit_depth == 12) {
        *full_lambda = av1lambda_mode_decision12_bit_sse[qp_index];
        *fast_lambda = av1lambda_mode_decision12_bit_sad[qp_index];
    }
    else {
        assert(bit_depth >= 8);
        assert(bit_depth <= 12);
    }

    //*full_lambda = 0; //-------------Nader
    //*fast_lambda = 0;
    *fast_chroma_lambda = *fast_lambda;
    *full_chroma_lambda = *full_lambda;

    // NM: To be done: tune lambda based on the picture type and layer.
}
const EbAv1LambdaAssignFunc av1_lambda_assignment_function_table[4] = {
    Av1lambdaAssign,
    Av1lambdaAssign,
    Av1lambdaAssign,
    Av1lambdaAssign,
};

void reset_mode_decision(
#if EIGHT_PEL_PREDICTIVE_ME
    SequenceControlSet    *sequence_control_set_ptr,
#endif
    ModeDecisionContext   *context_ptr,
    PictureControlSet     *picture_control_set_ptr,
    uint32_t                   segment_index)
{
    FrameHeader *frm_hdr = &picture_control_set_ptr->parent_pcs_ptr->frm_hdr;

    // QP
#if ADD_DELTA_QP_SUPPORT
    uint16_t picture_qp = picture_control_set_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx;
    context_ptr->qp = picture_qp;
    context_ptr->qp_index = context_ptr->qp;
#else
    context_ptr->qp = picture_control_set_ptr->picture_qp;
#endif
    // Asuming cb and cr offset to be the same for chroma QP in both slice and pps for lambda computation
    context_ptr->chroma_qp = (uint8_t)context_ptr->qp;
    context_ptr->qp_index = (uint8_t)frm_hdr->quantization_params.base_q_idx;
    (*av1_lambda_assignment_function_table[picture_control_set_ptr->parent_pcs_ptr->pred_structure])(
        &context_ptr->fast_lambda,
        &context_ptr->full_lambda,
        &context_ptr->fast_chroma_lambda,
        &context_ptr->full_chroma_lambda,
        (uint8_t)picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr->bit_depth,
        context_ptr->qp_index,
        context_ptr->hbd_mode_decision);
    // Reset MD rate Estimation table to initial values by copying from md_rate_estimation_array
    if (context_ptr->is_md_rate_estimation_ptr_owner) {
        context_ptr->is_md_rate_estimation_ptr_owner = EB_FALSE;
        EB_FREE_ARRAY(context_ptr->md_rate_estimation_ptr);
    }
    context_ptr->md_rate_estimation_ptr = picture_control_set_ptr->md_rate_estimation_array;
    uint32_t  candidateIndex;
    for (candidateIndex = 0; candidateIndex < MODE_DECISION_CANDIDATE_MAX_COUNT; ++candidateIndex)
        context_ptr->fast_candidate_ptr_array[candidateIndex]->md_rate_estimation_ptr = context_ptr->md_rate_estimation_ptr;

    // Reset CABAC Contexts
    context_ptr->coeff_est_entropy_coder_ptr = picture_control_set_ptr->coeff_est_entropy_coder_ptr;

    // Reset Neighbor Arrays at start of new Segment / Picture
    if (segment_index == 0) {
        reset_mode_decision_neighbor_arrays(picture_control_set_ptr);
#if !FIX_SETTINGS_RESET
    }
#endif
#if EIGHT_PEL_FIX
    (void)sequence_control_set_ptr;
#else
#if EIGHT_PEL_PREDICTIVE_ME
    picture_control_set_ptr->parent_pcs_ptr->frm_hdr.allow_high_precision_mv = picture_control_set_ptr->enc_mode == ENC_M0 &&
        (sequence_control_set_ptr->input_resolution == INPUT_SIZE_576p_RANGE_OR_LOWER) ? 1 : 0;
#else
#if EIGTH_PEL_MV
    picture_control_set_ptr->parent_pcs_ptr->allow_high_precision_mv = picture_control_set_ptr->enc_mode == ENC_M0 &&
        (picture_control_set_ptr->parent_pcs_ptr->is_pan || picture_control_set_ptr->parent_pcs_ptr->is_tilt) ? 1 : 0;
#endif
#endif
#endif
#if !MULTI_PASS_PD
    EbBool enable_wm;
    if (picture_control_set_ptr->parent_pcs_ptr->sc_content_detected)
        enable_wm = EB_FALSE;
    else
#if WARP_UPDATE
        enable_wm = (MR_MODE ||
        (picture_control_set_ptr->parent_pcs_ptr->enc_mode == ENC_M0 && picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag) ||
        (picture_control_set_ptr->parent_pcs_ptr->enc_mode <= ENC_M5 && picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index == 0)) ? EB_TRUE : EB_FALSE;
#else
        enable_wm = (picture_control_set_ptr->parent_pcs_ptr->enc_mode <= ENC_M5) || MR_MODE ? EB_TRUE : EB_FALSE;
#endif
#if !FIX_WM_SETTINGS
    enable_wm = picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index > 0 ? EB_FALSE : enable_wm;
#endif
    frm_hdr->allow_warped_motion = enable_wm
        && !(frm_hdr->frame_type == KEY_FRAME || frm_hdr->frame_type == INTRA_ONLY_FRAME)
        && !frm_hdr->error_resilient_mode;
    frm_hdr->is_motion_mode_switchable = frm_hdr->allow_warped_motion;
#if OBMC_FLAG
    // OBMC Level                                   Settings
    // 0                                            OFF
    // 1                                            OBMC @(MVP, PME and ME) + 16 NICs
    // 2                                            OBMC @(MVP, PME and ME) + Opt NICs
    // 3                                            OBMC @(MVP, PME ) + Opt NICs
    // 4                                            OBMC @(MVP, PME ) + Opt2 NICs
    if (sequence_control_set_ptr->static_config.enable_obmc) {
        if (picture_control_set_ptr->parent_pcs_ptr->enc_mode <= ENC_M0)
            picture_control_set_ptr->parent_pcs_ptr->pic_obmc_mode =
#if M0_OPT
            picture_control_set_ptr->slice_type != I_SLICE ? 2 : 0;
#else
            picture_control_set_ptr->parent_pcs_ptr->sc_content_detected == 0 && picture_control_set_ptr->slice_type != I_SLICE ? 2 : 0;
#endif
        else
            picture_control_set_ptr->parent_pcs_ptr->pic_obmc_mode = 0;

#if MR_MODE
        picture_control_set_ptr->parent_pcs_ptr->pic_obmc_mode =
            picture_control_set_ptr->parent_pcs_ptr->sc_content_detected == 0 && picture_control_set_ptr->slice_type != I_SLICE ? 1 : 0;
#endif
    }
    else
        picture_control_set_ptr->parent_pcs_ptr->pic_obmc_mode = 0;

    frm_hdr->is_motion_mode_switchable =
        frm_hdr->is_motion_mode_switchable || picture_control_set_ptr->parent_pcs_ptr->pic_obmc_mode;

#endif
#endif
#if FIX_SETTINGS_RESET
    }
#endif
    return;
}

/******************************************************
 * Mode Decision Configure LCU
 ******************************************************/
void mode_decision_configure_lcu(
    ModeDecisionContext   *context_ptr,
    PictureControlSet     *picture_control_set_ptr,
    uint8_t                    sb_qp){
    (void)picture_control_set_ptr;
    //Disable Lambda update per LCU
    context_ptr->qp = sb_qp;
    // Asuming cb and cr offset to be the same for chroma QP in both slice and pps for lambda computation

    context_ptr->chroma_qp = (uint8_t)context_ptr->qp;

    /* Note(CHKN) : when Qp modulation varies QP on a sub-LCU(CU) basis,  Lamda has to change based on Cu->QP , and then this code has to move inside the CU loop in MD */

    // Lambda Assignement
    context_ptr->qp_index = picture_control_set_ptr->parent_pcs_ptr->frm_hdr.delta_q_params.delta_q_present ? (uint8_t)quantizer_to_qindex[sb_qp] : (uint8_t)picture_control_set_ptr->parent_pcs_ptr->frm_hdr.quantization_params.base_q_idx;

    (*av1_lambda_assignment_function_table[picture_control_set_ptr->parent_pcs_ptr->pred_structure])(
        &context_ptr->fast_lambda,
        &context_ptr->full_lambda,
        &context_ptr->fast_chroma_lambda,
        &context_ptr->full_chroma_lambda,
        (uint8_t)picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr->bit_depth,
        context_ptr->qp_index,
        context_ptr->hbd_mode_decision);

    return;
}
