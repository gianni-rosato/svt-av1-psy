/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>

#include "EbUtility.h"
#include "EbModeDecisionProcess.h"
#include "EbLambdaRateTables.h"

/******************************************************
 * Mode Decision Context Constructor
 ******************************************************/
EbErrorType mode_decision_context_ctor(
    ModeDecisionContext  **context_dbl_ptr,
    EbColorFormat         color_format,
    EbFifo                *mode_decision_configuration_input_fifo_ptr,
    EbFifo                *mode_decision_output_fifo_ptr){
    uint32_t bufferIndex;
    uint32_t candidateIndex;
    EbErrorType return_error = EB_ErrorNone;

    (void)color_format;

    ModeDecisionContext *context_ptr;
    EB_MALLOC(ModeDecisionContext*, context_ptr, sizeof(ModeDecisionContext), EB_N_PTR);
    *context_dbl_ptr = context_ptr;

    // Input/Output System Resource Manager FIFOs
    context_ptr->mode_decision_configuration_input_fifo_ptr = mode_decision_configuration_input_fifo_ptr;
    context_ptr->mode_decision_output_fifo_ptr = mode_decision_output_fifo_ptr;

    // Trasform Scratch Memory
    EB_MALLOC(int16_t*, context_ptr->transform_inner_array_ptr, 3120, EB_N_PTR); //refer to EbInvTransform_SSE2.as. case 32x32

    // MD rate Estimation tables
    EB_MALLOC(MdRateEstimationContext*, context_ptr->md_rate_estimation_ptr, sizeof(MdRateEstimationContext), EB_N_PTR);

    // Fast Candidate Array
    EB_MALLOC(ModeDecisionCandidate*, context_ptr->fast_candidate_array, sizeof(ModeDecisionCandidate) * MODE_DECISION_CANDIDATE_MAX_COUNT, EB_N_PTR);

    EB_MALLOC(ModeDecisionCandidate**, context_ptr->fast_candidate_ptr_array, sizeof(ModeDecisionCandidate*) * MODE_DECISION_CANDIDATE_MAX_COUNT, EB_N_PTR);

    for (candidateIndex = 0; candidateIndex < MODE_DECISION_CANDIDATE_MAX_COUNT; ++candidateIndex) {
        context_ptr->fast_candidate_ptr_array[candidateIndex] = &context_ptr->fast_candidate_array[candidateIndex];
        context_ptr->fast_candidate_ptr_array[candidateIndex]->md_rate_estimation_ptr = context_ptr->md_rate_estimation_ptr;
    }

    // Transform and Quantization Buffers
    EB_MALLOC(EbTransQuantBuffers*, context_ptr->trans_quant_buffers_ptr, sizeof(EbTransQuantBuffers), EB_N_PTR);

    return_error = eb_trans_quant_buffers_ctor(
        context_ptr->trans_quant_buffers_ptr);

    if (return_error == EB_ErrorInsufficientResources)
        return EB_ErrorInsufficientResources;
    // Cost Arrays
    // Hsan: MAX_NFL + 1 scratch buffer for intra + 1 scratch buffer for inter
    EB_MALLOC(uint64_t*, context_ptr->fast_cost_array, sizeof(uint64_t) * (MAX_NFL + 1 + 1), EB_N_PTR);
    EB_MALLOC(uint64_t*, context_ptr->full_cost_array, sizeof(uint64_t) * (MAX_NFL + 1 + 1), EB_N_PTR);
    EB_MALLOC(uint64_t*, context_ptr->full_cost_skip_ptr, sizeof(uint64_t) * (MAX_NFL + 1 + 1), EB_N_PTR);
    EB_MALLOC(uint64_t*, context_ptr->full_cost_merge_ptr, sizeof(uint64_t) * (MAX_NFL + 1 + 1), EB_N_PTR);
    // Candidate Buffers
    EB_MALLOC(ModeDecisionCandidateBuffer**, context_ptr->candidate_buffer_ptr_array, sizeof(ModeDecisionCandidateBuffer*) * (MAX_NFL + 1 + 1), EB_N_PTR);

    for (bufferIndex = 0; bufferIndex < (MAX_NFL + 1 + 1); ++bufferIndex) {
        return_error = mode_decision_candidate_buffer_ctor(
            &(context_ptr->candidate_buffer_ptr_array[bufferIndex]),
            &(context_ptr->fast_cost_array[bufferIndex]),
            &(context_ptr->full_cost_array[bufferIndex]),
            &(context_ptr->full_cost_skip_ptr[bufferIndex]),
            &(context_ptr->full_cost_merge_ptr[bufferIndex])
        );
        if (return_error == EB_ErrorInsufficientResources)
            return EB_ErrorInsufficientResources;
    }
#if !UNPACK_REF_POST_EP
    // Inter Prediction Context
    return_error = inter_prediction_context_ctor(
        &context_ptr->inter_prediction_context,
        color_format,
        SB_STRIDE_Y,
        SB_STRIDE_Y);
    if (return_error == EB_ErrorInsufficientResources)
        return EB_ErrorInsufficientResources;
#endif
    // Intra Reference Samples
    return_error = intra_reference_samples_ctor(&context_ptr->intra_ref_ptr);
    if (return_error == EB_ErrorInsufficientResources)
        return EB_ErrorInsufficientResources;
    uint32_t codedLeafIndex, tu_index;

    for (codedLeafIndex = 0; codedLeafIndex < BLOCK_MAX_COUNT_SB_128; ++codedLeafIndex) {
        for (tu_index = 0; tu_index < TRANSFORM_UNIT_MAX_COUNT; ++tu_index)
            context_ptr->md_cu_arr_nsq[codedLeafIndex].transform_unit_array[tu_index].tu_index = tu_index;
        const BlockGeom * blk_geom = get_blk_geom_mds(codedLeafIndex);
        UNUSED(blk_geom);
        EB_MALLOC(MacroBlockD*, context_ptr->md_cu_arr_nsq[codedLeafIndex].av1xd, sizeof(MacroBlockD), EB_N_PTR);
        EB_MALLOC(uint8_t*, context_ptr->md_cu_arr_nsq[codedLeafIndex].neigh_left_recon[0], 128, EB_N_PTR);
        EB_MALLOC(uint8_t*, context_ptr->md_cu_arr_nsq[codedLeafIndex].neigh_left_recon[1], 128, EB_N_PTR);
        EB_MALLOC(uint8_t*, context_ptr->md_cu_arr_nsq[codedLeafIndex].neigh_left_recon[2], 128, EB_N_PTR);

        EB_MALLOC(uint8_t*, context_ptr->md_cu_arr_nsq[codedLeafIndex].neigh_top_recon[0], 128, EB_N_PTR);
        EB_MALLOC(uint8_t*, context_ptr->md_cu_arr_nsq[codedLeafIndex].neigh_top_recon[1], 128, EB_N_PTR);
        EB_MALLOC(uint8_t*, context_ptr->md_cu_arr_nsq[codedLeafIndex].neigh_top_recon[2], 128, EB_N_PTR);

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

            return_error = eb_picture_buffer_desc_ctor(
                (EbPtr*)&context_ptr->md_cu_arr_nsq[codedLeafIndex].coeff_tmp,
                (EbPtr)&initData);

            if (return_error == EB_ErrorInsufficientResources)
                return EB_ErrorInsufficientResources;
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

            return_error = eb_picture_buffer_desc_ctor(
                (EbPtr*)&context_ptr->md_cu_arr_nsq[codedLeafIndex].recon_tmp,
                (EbPtr)&initData);

            if (return_error == EB_ErrorInsufficientResources)
                return EB_ErrorInsufficientResources;
        }
#endif
    }
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

        neighbor_array_unit_reset(picture_control_set_ptr->md_luma_recon_neighbor_array[depth]);
#if ATB_MD
        neighbor_array_unit_reset(picture_control_set_ptr->md_tx_depth_1_luma_recon_neighbor_array[depth]);
#endif
        neighbor_array_unit_reset(picture_control_set_ptr->md_cb_recon_neighbor_array[depth]);
        neighbor_array_unit_reset(picture_control_set_ptr->md_cr_recon_neighbor_array[depth]);
#if !REMOVE_SKIP_COEFF_NEIGHBOR_ARRAY
        neighbor_array_unit_reset(picture_control_set_ptr->md_skip_coeff_neighbor_array[depth]);
#endif
        neighbor_array_unit_reset(picture_control_set_ptr->md_luma_dc_sign_level_coeff_neighbor_array[depth]);
#if ATB_DC_CONTEXT_SUPPORT_2
        neighbor_array_unit_reset(picture_control_set_ptr->md_tx_depth_1_luma_dc_sign_level_coeff_neighbor_array[depth]);
#endif
        neighbor_array_unit_reset(picture_control_set_ptr->md_cb_dc_sign_level_coeff_neighbor_array[depth]);
        neighbor_array_unit_reset(picture_control_set_ptr->md_cr_dc_sign_level_coeff_neighbor_array[depth]);
#if ATB_RATE
        neighbor_array_unit_reset(picture_control_set_ptr->md_txfm_context_array[depth]);
#endif
        neighbor_array_unit_reset(picture_control_set_ptr->md_inter_pred_dir_neighbor_array[depth]);
        neighbor_array_unit_reset(picture_control_set_ptr->md_ref_frame_type_neighbor_array[depth]);

        neighbor_array_unit_reset32(picture_control_set_ptr->md_interpolation_type_neighbor_array[depth]);
    }

    return;
}

#if !MEMORY_FOOTPRINT_OPT
void ResetMdRefinmentNeighborArrays(PictureControlSet *picture_control_set_ptr)
{
    neighbor_array_unit_reset(picture_control_set_ptr->md_refinement_mode_type_neighbor_array);
    neighbor_array_unit_reset(picture_control_set_ptr->md_refinement_luma_recon_neighbor_array);
    return;
}
#endif

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
    uint16_t                     qp_index)

{
    if (bit_depth == 8) {
        *full_lambda = av1_lambda_mode_decision8_bit_sse[qp_index];
        *fast_lambda = av1_lambda_mode_decision8_bit_sad[qp_index];
    }
    else if (bit_depth == 10) {
        *full_lambda = av1lambda_mode_decision10_bit_sse[qp_index];
        *fast_lambda = av1lambda_mode_decision10_bit_sad[qp_index];
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
    ModeDecisionContext   *context_ptr,
    PictureControlSet     *picture_control_set_ptr,
    SequenceControlSet    *sequence_control_set_ptr,
    uint32_t                   segment_index)
{
    EB_SLICE                     slice_type;
#if !MEMORY_FOOTPRINT_OPT
    uint32_t                       lcuRowIndex;
#endif
    MdRateEstimationContext   *md_rate_estimation_array;

    // QP
#if ADD_DELTA_QP_SUPPORT
    uint16_t picture_qp = picture_control_set_ptr->parent_pcs_ptr->base_qindex;
    context_ptr->qp = picture_qp;
    context_ptr->qp_index = context_ptr->qp;
#else
    context_ptr->qp = picture_control_set_ptr->picture_qp;
#endif
    // Asuming cb and cr offset to be the same for chroma QP in both slice and pps for lambda computation
    context_ptr->chroma_qp = context_ptr->qp;
    context_ptr->qp_index = (uint8_t)picture_control_set_ptr->parent_pcs_ptr->base_qindex;
    (*av1_lambda_assignment_function_table[picture_control_set_ptr->parent_pcs_ptr->pred_structure])(
        &context_ptr->fast_lambda,
        &context_ptr->full_lambda,
        &context_ptr->fast_chroma_lambda,
        &context_ptr->full_chroma_lambda,
        (uint8_t)picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr->bit_depth,
        context_ptr->qp_index);
    // Slice Type
    slice_type =
        (picture_control_set_ptr->parent_pcs_ptr->idr_flag == EB_TRUE) ? I_SLICE :
        picture_control_set_ptr->slice_type;

    // Increment the MD Rate Estimation array pointer to point to the right address based on the QP and slice type

    /* Note(CHKN) : Rate estimation will use FrameQP even when Qp modulation is ON */

    md_rate_estimation_array = (MdRateEstimationContext*)sequence_control_set_ptr->encode_context_ptr->md_rate_estimation_array;
#if ADD_DELTA_QP_SUPPORT
    md_rate_estimation_array += slice_type * TOTAL_NUMBER_OF_QP_VALUES + picture_control_set_ptr->picture_qp;
#else
    md_rate_estimation_array += slice_type * TOTAL_NUMBER_OF_QP_VALUES + context_ptr->qp;
#endif

    // Reset MD rate Estimation table to initial values by copying from md_rate_estimation_array

    context_ptr->md_rate_estimation_ptr = md_rate_estimation_array;
    uint32_t  candidateIndex;
    for (candidateIndex = 0; candidateIndex < MODE_DECISION_CANDIDATE_MAX_COUNT; ++candidateIndex)
        context_ptr->fast_candidate_ptr_array[candidateIndex]->md_rate_estimation_ptr = md_rate_estimation_array;
#if !OPT_LOSSLESS_0
    // TMVP Map Writer pointer
    if (picture_control_set_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE)
        context_ptr->reference_object_write_ptr = (EbReferenceObject*)picture_control_set_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr;
    else
        context_ptr->reference_object_write_ptr = (EbReferenceObject*)EB_NULL;
#endif
    // Reset CABAC Contexts
    context_ptr->coeff_est_entropy_coder_ptr = picture_control_set_ptr->coeff_est_entropy_coder_ptr;

    // Reset Neighbor Arrays at start of new Segment / Picture
    if (segment_index == 0) {
        reset_mode_decision_neighbor_arrays(picture_control_set_ptr);
#if !MEMORY_FOOTPRINT_OPT
        ResetMdRefinmentNeighborArrays(picture_control_set_ptr);
        for (lcuRowIndex = 0; lcuRowIndex < ((sequence_control_set_ptr->luma_height + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64); lcuRowIndex++) {
            picture_control_set_ptr->enc_prev_coded_qp[lcuRowIndex] = (uint8_t)picture_control_set_ptr->picture_qp;
            picture_control_set_ptr->enc_prev_quant_group_coded_qp[lcuRowIndex] = (uint8_t)picture_control_set_ptr->picture_qp;
        }
#endif
    }

#if EIGTH_PEL_MV
    picture_control_set_ptr->parent_pcs_ptr->allow_high_precision_mv = picture_control_set_ptr->enc_mode == ENC_M0 &&
        (picture_control_set_ptr->parent_pcs_ptr->is_pan || picture_control_set_ptr->parent_pcs_ptr->is_tilt) ? 1 : 0;
#endif

#if ENABLE_WARPED_MV
#if NEW_PRESETS
    EbBool enable_wm = (picture_control_set_ptr->parent_pcs_ptr->enc_mode <= ENC_M5) || MR_MODE ? EB_TRUE : EB_FALSE;
#else
    EbBool enable_wm = (picture_control_set_ptr->parent_pcs_ptr->enc_mode == ENC_M0) || MR_MODE ? EB_TRUE : EB_FALSE;
#endif
    enable_wm = picture_control_set_ptr->parent_pcs_ptr->temporal_layer_index > 0 ? EB_FALSE : enable_wm;
    picture_control_set_ptr->parent_pcs_ptr->allow_warped_motion = enable_wm
#else
    picture_control_set_ptr->parent_pcs_ptr->allow_warped_motion = sequence_control_set_ptr->static_config.enable_warped_motion
#endif
        && !(picture_control_set_ptr->parent_pcs_ptr->av1_frame_type == KEY_FRAME || picture_control_set_ptr->parent_pcs_ptr->av1_frame_type == INTRA_ONLY_FRAME)
        && !picture_control_set_ptr->parent_pcs_ptr->error_resilient_mode;
    picture_control_set_ptr->parent_pcs_ptr->switchable_motion_mode = picture_control_set_ptr->parent_pcs_ptr->allow_warped_motion;

    return;
}

/******************************************************
 * Mode Decision Configure LCU
 ******************************************************/
void mode_decision_configure_lcu(
    ModeDecisionContext   *context_ptr,
    LargestCodingUnit     *sb_ptr,
    PictureControlSet     *picture_control_set_ptr,
    SequenceControlSet    *sequence_control_set_ptr,
    uint8_t                    picture_qp,
    uint8_t                    sb_qp){
    (void)picture_control_set_ptr;
    //Disable Lambda update per LCU

    //RC is off
    if (sequence_control_set_ptr->static_config.rate_control_mode == 0 && sequence_control_set_ptr->static_config.improve_sharpness == 0) {
        context_ptr->qp = (uint8_t)picture_qp;
        sb_ptr->qp = (uint8_t)context_ptr->qp;
    }
    //RC is on
    else
        context_ptr->qp = (uint8_t)sb_qp;
    // Asuming cb and cr offset to be the same for chroma QP in both slice and pps for lambda computation

    context_ptr->chroma_qp = context_ptr->qp;

    /* Note(CHKN) : when Qp modulation varies QP on a sub-LCU(CU) basis,  Lamda has to change based on Cu->QP , and then this code has to move inside the CU loop in MD */

    // Lambda Assignement
    context_ptr->qp_index = (uint8_t)picture_control_set_ptr->parent_pcs_ptr->base_qindex;

    (*av1_lambda_assignment_function_table[picture_control_set_ptr->parent_pcs_ptr->pred_structure])(
        &context_ptr->fast_lambda,
        &context_ptr->full_lambda,
        &context_ptr->fast_chroma_lambda,
        &context_ptr->full_chroma_lambda,
        (uint8_t)picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr->bit_depth,
        context_ptr->qp_index);

    return;
}
