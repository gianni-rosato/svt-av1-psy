/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>

#include "EbDefinitions.h"
#include "EbSequenceControlSet.h"

/**************************************************************************************************
    General notes on how Sequence Control Sets (SCS) are used.

    SequenceControlSetInstance
        is the master copy that interacts with the API in real-time.  When a
        change happens, the changeFlag is signaled so that appropriate action can
        be taken.  There is one scsInstance per stream/encode instance.  The scsInstance
        owns the encodeContext

    encodeContext
        has context type variables (i.e. non-config) that keep track of global parameters.

    SequenceControlSets
        general SCSs are controled by a system resource manager.  They are kept completely
        separate from the instances.  In general there is one active SCS at a time.  When the
        changeFlag is signaled, the old active SCS is no longer used for new input pictures.
        A fresh copy of the scsInstance is made to a new SCS, which becomes the active SCS.  The
        old SCS will eventually be released back into the SCS pool when its current pictures are
        finished encoding.

    Motivations
        The whole reason for this structure is due to the nature of the pipeline.  We have to
        take great care not to have pipeline mismanagement.  Once an object enters use in the
        pipeline, it cannot be changed on the fly or you will have pipeline coherency problems.
 ***************************************************************************************************/
EbErrorType eb_sequence_control_set_ctor(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    EbSequenceControlSetInitData_t *scsInitData = (EbSequenceControlSetInitData_t*)object_init_data_ptr;
    uint32_t segment_index;
    SequenceControlSet_t *sequence_control_set_ptr;
    EB_MALLOC(SequenceControlSet_t*, sequence_control_set_ptr, sizeof(SequenceControlSet_t), EB_N_PTR);

    *object_dbl_ptr = (EbPtr)sequence_control_set_ptr;

    sequence_control_set_ptr->static_config.sb_sz = 64;
    sequence_control_set_ptr->static_config.partition_depth = 4;
    sequence_control_set_ptr->static_config.qp = 32;

    // Segments
    for (segment_index = 0; segment_index < MAX_TEMPORAL_LAYERS; ++segment_index) {
        sequence_control_set_ptr->me_segment_column_count_array[segment_index] = 1;
        sequence_control_set_ptr->me_segment_row_count_array[segment_index] = 1;
        sequence_control_set_ptr->enc_dec_segment_col_count_array[segment_index] = 1;
        sequence_control_set_ptr->enc_dec_segment_row_count_array[segment_index] = 1;
    }

    // Encode Context
    if (scsInitData != EB_NULL) {
        sequence_control_set_ptr->encode_context_ptr = scsInitData->encode_context_ptr;
    }
    else {
        sequence_control_set_ptr->encode_context_ptr = (EncodeContext_t *)EB_NULL;
    }

    sequence_control_set_ptr->conformance_window_flag = 0;

    // Profile & ID
    sequence_control_set_ptr->sps_id = 0;
    sequence_control_set_ptr->vps_id = 0;
    sequence_control_set_ptr->profile_space = 0;
    sequence_control_set_ptr->profile_idc = 0;
    sequence_control_set_ptr->level_idc = 0;
    sequence_control_set_ptr->tier_idc = 0;
    sequence_control_set_ptr->chroma_format_idc = 1; // EB_YUV420
    sequence_control_set_ptr->max_temporal_layers = 1;

    sequence_control_set_ptr->bits_for_picture_order_count = 16;
    sequence_control_set_ptr->subsampling_y = 1;
    sequence_control_set_ptr->subsampling_x = 1;

    // Picture Dimensions
    sequence_control_set_ptr->luma_width = 0;
    sequence_control_set_ptr->luma_height = 0;

    sequence_control_set_ptr->chroma_width = 0;
    sequence_control_set_ptr->chroma_height = 0;
    sequence_control_set_ptr->frame_rate = 0;
    sequence_control_set_ptr->encoder_bit_depth = 8;

    // Bitdepth
    sequence_control_set_ptr->input_bitdepth = EB_8BIT;
    sequence_control_set_ptr->output_bitdepth = EB_8BIT;

    // GOP Structure
    sequence_control_set_ptr->max_ref_count = 1;
    sequence_control_set_ptr->intra_period_length = 0;
    sequence_control_set_ptr->intra_refresh_type = 0;

    // LCU
    sequence_control_set_ptr->sb_sz = 64;
    sequence_control_set_ptr->max_sb_depth = 3;

    // CU
    sequence_control_set_ptr->max_cu_size = 64;
    sequence_control_set_ptr->min_cu_size = 8;
    sequence_control_set_ptr->max_intra_size = 32;
    sequence_control_set_ptr->min_intra_size = 8;
    sequence_control_set_ptr->intra4x4_flag = EB_TRUE;

    sequence_control_set_ptr->general_progressive_source_flag = 1;
    sequence_control_set_ptr->general_interlaced_source_flag = 0;
    sequence_control_set_ptr->general_frame_only_constraint_flag = 0;

    // temporal mvp enable flag
    sequence_control_set_ptr->enable_tmvp_sps = 1;

    // Rate Control
    sequence_control_set_ptr->rate_control_mode = 0;
    sequence_control_set_ptr->target_bitrate = 0x1000;
    sequence_control_set_ptr->available_bandwidth = 0x1000;

    // Quantization
    sequence_control_set_ptr->qp = 20;

    // mv merge
    sequence_control_set_ptr->mv_merge_total_count = 5;

    // Initialize SB params
    sb_params_ctor(
        sequence_control_set_ptr);

    sequence_control_set_ptr->num_bits_width = 16;
    sequence_control_set_ptr->num_bits_height = 16;
    sequence_control_set_ptr->frame_id_numbers_present_flag = 0;
    //    cm->large_scale_tile ? 0 : cm->error_resilient_mode;
    sequence_control_set_ptr->frame_id_length = FRAME_ID_LENGTH;
    sequence_control_set_ptr->delta_frame_id_length = DELTA_FRAME_ID_LENGTH;

    if (scsInitData && scsInitData->sb_size == 128)
    {
        sequence_control_set_ptr->sb_size = BLOCK_128X128;
        sequence_control_set_ptr->sb_size_pix = 128;
        sequence_control_set_ptr->max_block_cnt = 4421;

        sequence_control_set_ptr->mib_size = 32;        // Size of the superblock in units of MI blocks
        sequence_control_set_ptr->mib_size_log2 = 5;
    }
    else
    {
        sequence_control_set_ptr->sb_size = BLOCK_64X64;
        sequence_control_set_ptr->sb_size_pix = 64;
        sequence_control_set_ptr->max_block_cnt = 1101;

        sequence_control_set_ptr->mib_size = 16;        // Size of the superblock in units of MI blocks
        sequence_control_set_ptr->mib_size_log2 = 4;
    }



    sequence_control_set_ptr->enable_dual_filter = 1;
    // 0 - disable dual interpolation filter
    // 1 - enable vertical and horiz filter selection
    sequence_control_set_ptr->enable_order_hint = 1;
    // 0 - disable order hint, and related tools:
    // jnt_comp, ref_frame_mvs, frame_sign_bias
    // if 0, enable_jnt_comp must be set zs 0.
    sequence_control_set_ptr->enable_jnt_comp = 0;

    sequence_control_set_ptr->order_hint_bits_minus1 = sequence_control_set_ptr->enable_order_hint ? 6 : -1;

    sequence_control_set_ptr->force_screen_content_tools = 0;
    // 0 - force off
    // 1 - force on
    // 2 - adaptive
    sequence_control_set_ptr->force_integer_mv = 2;  // 0 - Not to force. MV can be in 1/4 or 1/8
    // 1 - force to integer
    // 2 - adaptive

    sequence_control_set_ptr->enable_filter_intra = 0;
    sequence_control_set_ptr->enable_intra_edge_filter = 0;

    sequence_control_set_ptr->enable_interintra_compound = 0;
    sequence_control_set_ptr->enable_masked_compound = 0;

    sequence_control_set_ptr->enable_ref_frame_mvs = 1;
    sequence_control_set_ptr->enable_superres = 0;
#if NO_ENCDEC || SHUT_FILTERING
    sequence_control_set_ptr->enable_cdef = 0;

    sequence_control_set_ptr->enable_restoration = 0;
#else
    sequence_control_set_ptr->enable_cdef = 1;

    sequence_control_set_ptr->film_grain_params_present = 0;
    sequence_control_set_ptr->film_grain_denoise_strength = 0;

    sequence_control_set_ptr->enable_restoration = 1;
#endif
    sequence_control_set_ptr->reduced_still_picture_hdr = 0;
    sequence_control_set_ptr->still_picture = 0;
    sequence_control_set_ptr->timing_info_present = 0;
    sequence_control_set_ptr->operating_points_decoder_model_cnt = 0;
#if AV1_UPGRADE
    sequence_control_set_ptr->decoder_model_info_present_flag = 0;
    sequence_control_set_ptr->display_model_info_present_flag = 0;
#endif
    for (int32_t i = 0; i < MAX_NUM_OPERATING_POINTS; i++) {
        sequence_control_set_ptr->operating_point_idc[i] = 0;
        sequence_control_set_ptr->level[i].major = 0;
        sequence_control_set_ptr->level[i].minor = 0;
        sequence_control_set_ptr->tier[i] = 0;
    }
    sequence_control_set_ptr->monochrome = 0;
    sequence_control_set_ptr->film_grain_params_present = 0;
    sequence_control_set_ptr->film_grain_random_seed = 7391;
    return EB_ErrorNone;
}


/************************************************
 * Sequence Control Set Copy
 ************************************************/
EbErrorType copy_sequence_control_set(
    SequenceControlSet_t *dst,
    SequenceControlSet_t *src)
{
    uint32_t  writeCount = 0;

    dst->static_config = src->static_config;                            writeCount += sizeof(EbSvtAv1EncConfiguration);
    dst->encode_context_ptr = src->encode_context_ptr;                        writeCount += sizeof(EncodeContext_t*);
    dst->sps_id = src->sps_id;                                   writeCount += sizeof(uint32_t);
    dst->vps_id = src->vps_id;                                   writeCount += sizeof(uint32_t);
    dst->profile_space = src->profile_space;                            writeCount += sizeof(uint32_t);
    dst->profile_idc = src->profile_idc;                              writeCount += sizeof(uint32_t);
    dst->level_idc = src->level_idc;                                writeCount += sizeof(uint32_t);
    dst->tier_idc = src->tier_idc;                                 writeCount += sizeof(uint32_t);
    dst->chroma_format_idc = src->chroma_format_idc;                         writeCount += sizeof(uint32_t);
    dst->max_temporal_layers = src->max_temporal_layers;                       writeCount += sizeof(uint32_t);
    dst->bits_for_picture_order_count = src->bits_for_picture_order_count;                writeCount += sizeof(uint32_t);
    dst->max_input_luma_width = src->max_input_luma_width;                       writeCount += sizeof(uint32_t);
    dst->max_input_luma_height = src->max_input_luma_height;                      writeCount += sizeof(uint32_t);
    dst->max_input_chroma_height = src->max_input_chroma_height;                    writeCount += sizeof(uint32_t);
    dst->max_input_chroma_width = src->max_input_chroma_width;                     writeCount += sizeof(uint32_t);
    dst->max_input_pad_right = src->max_input_pad_right;                        writeCount += sizeof(uint32_t);
    dst->max_input_pad_bottom = src->max_input_pad_bottom;                       writeCount += sizeof(uint32_t);
    dst->luma_width = src->luma_width;                               writeCount += sizeof(uint32_t);
    dst->luma_height = src->luma_height;                              writeCount += sizeof(uint32_t);
    dst->chroma_width = src->chroma_width;                             writeCount += sizeof(uint32_t);
    dst->chroma_height = src->chroma_height;                            writeCount += sizeof(uint32_t);
    dst->pad_right = src->pad_right;                                writeCount += sizeof(uint32_t);
    dst->pad_bottom = src->pad_bottom;                               writeCount += sizeof(uint32_t);
    dst->cropping_left_offset = src->cropping_left_offset;                      writeCount += sizeof(int32_t);
    dst->cropping_right_offset = src->cropping_right_offset;                     writeCount += sizeof(int32_t);
    dst->cropping_top_offset = src->cropping_top_offset;                       writeCount += sizeof(int32_t);
    dst->cropping_bottom_offset = src->cropping_bottom_offset;                    writeCount += sizeof(int32_t);
    dst->conformance_window_flag = src->conformance_window_flag;                   writeCount += sizeof(uint32_t);
    dst->frame_rate = src->frame_rate;                               writeCount += sizeof(uint32_t);
    dst->input_bitdepth = src->input_bitdepth;                           writeCount += sizeof(EB_BITDEPTH);
    dst->output_bitdepth = src->output_bitdepth;                          writeCount += sizeof(EB_BITDEPTH);
    dst->pred_struct_ptr = src->pred_struct_ptr;                           writeCount += sizeof(PredictionStructure_t*);
    dst->intra_period_length = src->intra_period_length;                       writeCount += sizeof(int32_t);
    dst->intra_refresh_type = src->intra_refresh_type;                        writeCount += sizeof(uint32_t);
    dst->max_ref_count = src->max_ref_count;                             writeCount += sizeof(uint32_t);
    dst->sb_sz = src->sb_sz;                                 writeCount += sizeof(uint32_t);
    dst->max_sb_depth = src->max_sb_depth;                             writeCount += sizeof(uint32_t);
    dst->max_cu_size = src->max_cu_size;                               writeCount += sizeof(uint32_t);
    dst->min_cu_size = src->min_cu_size;                               writeCount += sizeof(uint32_t);
    dst->max_intra_size = src->max_intra_size;                            writeCount += sizeof(uint32_t);
    dst->min_intra_size = src->min_intra_size;                            writeCount += sizeof(uint32_t);
    dst->intra4x4_flag = src->intra4x4_flag;                            writeCount += sizeof(EbBool);
    dst->general_progressive_source_flag = src->general_progressive_source_flag;        writeCount += sizeof(uint32_t);
    dst->general_interlaced_source_flag = src->general_interlaced_source_flag;         writeCount += sizeof(uint32_t);
    dst->general_frame_only_constraint_flag = src->general_frame_only_constraint_flag;      writeCount += sizeof(uint32_t);
    dst->target_bitrate = src->target_bitrate;                           writeCount += sizeof(uint32_t);
    dst->available_bandwidth = src->available_bandwidth;                      writeCount += sizeof(uint32_t);
    dst->qp = src->qp;                                      writeCount += sizeof(uint32_t);
    dst->enable_tmvp_sps = src->enable_tmvp_sps;                           writeCount += sizeof(uint32_t);
    dst->mv_merge_total_count = src->mv_merge_total_count;                       writeCount += sizeof(uint32_t);
    dst->film_grain_denoise_strength = src->film_grain_denoise_strength;          writeCount += sizeof(int32_t);
    dst->film_grain_params_present = src->film_grain_params_present;              writeCount += sizeof(int32_t);
    dst->film_grain_params_present = src->film_grain_params_present;              writeCount += sizeof(int32_t);
    dst->picture_control_set_pool_init_count = src->picture_control_set_pool_init_count;            writeCount += sizeof(int32_t);
    dst->picture_control_set_pool_init_count_child = src->picture_control_set_pool_init_count_child; writeCount += sizeof(int32_t);
    dst->pa_reference_picture_buffer_init_count = src->pa_reference_picture_buffer_init_count; writeCount += sizeof(int32_t);
    dst->reference_picture_buffer_init_count = src->reference_picture_buffer_init_count; writeCount += sizeof(int32_t);
    dst->input_buffer_fifo_init_count = src->input_buffer_fifo_init_count; writeCount += sizeof(int32_t);
    dst->output_stream_buffer_fifo_init_count = src->output_stream_buffer_fifo_init_count; writeCount += sizeof(int32_t);
    dst->output_recon_buffer_fifo_init_count = src->output_recon_buffer_fifo_init_count; writeCount += sizeof(int32_t);
    dst->resource_coordination_fifo_init_count = src->resource_coordination_fifo_init_count; writeCount += sizeof(int32_t);
    dst->picture_analysis_fifo_init_count = src->picture_analysis_fifo_init_count; writeCount += sizeof(int32_t);
    dst->picture_decision_fifo_init_count = src->picture_decision_fifo_init_count; writeCount += sizeof(int32_t);
    dst->motion_estimation_fifo_init_count = src->motion_estimation_fifo_init_count; writeCount += sizeof(int32_t);
    dst->initial_rate_control_fifo_init_count = src->initial_rate_control_fifo_init_count; writeCount += sizeof(int32_t);
    dst->picture_demux_fifo_init_count = src->picture_demux_fifo_init_count; writeCount += sizeof(int32_t);
    dst->rate_control_tasks_fifo_init_count = src->rate_control_tasks_fifo_init_count; writeCount += sizeof(int32_t);
    dst->rate_control_fifo_init_count = src->rate_control_fifo_init_count; writeCount += sizeof(int32_t);
    dst->mode_decision_configuration_fifo_init_count = src->mode_decision_configuration_fifo_init_count; writeCount += sizeof(int32_t);
    dst->enc_dec_fifo_init_count = src->enc_dec_fifo_init_count; writeCount += sizeof(int32_t);
    dst->entropy_coding_fifo_init_count = src->entropy_coding_fifo_init_count; writeCount += sizeof(int32_t);
    dst->picture_analysis_process_init_count = src->picture_analysis_process_init_count; writeCount += sizeof(int32_t);
    dst->motion_estimation_process_init_count = src->motion_estimation_process_init_count; writeCount += sizeof(int32_t);
    dst->source_based_operations_process_init_count = src->source_based_operations_process_init_count; writeCount += sizeof(int32_t);
    dst->mode_decision_configuration_process_init_count = src->mode_decision_configuration_process_init_count; writeCount += sizeof(int32_t);
    dst->enc_dec_process_init_count = src->enc_dec_process_init_count; writeCount += sizeof(int32_t);
    dst->entropy_coding_process_init_count = src->entropy_coding_process_init_count; writeCount += sizeof(int32_t);
    dst->total_process_init_count = src->total_process_init_count; writeCount += sizeof(int32_t);
    dst->left_padding = src->left_padding; writeCount += sizeof(int16_t);
    dst->right_padding = src->right_padding; writeCount += sizeof(int16_t);
    dst->top_padding = src->top_padding; writeCount += sizeof(int16_t);
    dst->bot_padding = src->bot_padding; writeCount += sizeof(int16_t);

    for (uint8_t i = 0; i< MAX_HIERARCHICAL_LEVEL; i++) {
        dst->me_segment_column_count_array[i] = src->me_segment_column_count_array[i];
        dst->me_segment_row_count_array[i] = src->me_segment_row_count_array[i];
        dst->enc_dec_segment_col_count_array[i] = src->enc_dec_segment_col_count_array[i];
        dst->enc_dec_segment_row_count_array[i] = src->enc_dec_segment_row_count_array[i];
    }

#if CDEF_M
    dst->cdef_segment_column_count = src->cdef_segment_column_count;
    dst->cdef_segment_row_count = src->cdef_segment_row_count;
#endif
#if REST_M
    dst->rest_segment_column_count = src->rest_segment_column_count;
    dst->rest_segment_row_count = src->rest_segment_row_count;
#endif
    return EB_ErrorNone;
}

extern EbErrorType derive_input_resolution(
    SequenceControlSet_t *sequenceControlSetPtr,
    uint32_t                  inputSize) {
    EbErrorType return_error = EB_ErrorNone;

    sequenceControlSetPtr->input_resolution = (inputSize < INPUT_SIZE_1080i_TH) ?
        INPUT_SIZE_576p_RANGE_OR_LOWER :
        (inputSize < INPUT_SIZE_1080p_TH) ?
        INPUT_SIZE_1080i_RANGE :
        (inputSize < INPUT_SIZE_4K_TH) ?
        INPUT_SIZE_1080p_RANGE :
        INPUT_SIZE_4K_RANGE;

    return return_error;
}

EbErrorType eb_sequence_control_set_instance_ctor(
    EbSequenceControlSetInstance_t **object_dbl_ptr)
{
    EbSequenceControlSetInitData_t scsInitData;
    EbErrorType return_error = EB_ErrorNone;
    EB_MALLOC(EbSequenceControlSetInstance_t*, *object_dbl_ptr, sizeof(EbSequenceControlSetInstance_t), EB_N_PTR);

    scsInitData.sb_size = 64;

    return_error = encode_context_ctor(
        (void **) &(*object_dbl_ptr)->encode_context_ptr,
        EB_NULL);
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }
    scsInitData.encode_context_ptr = (*object_dbl_ptr)->encode_context_ptr;

    scsInitData.sb_size = 64;

    return_error = eb_sequence_control_set_ctor(
        (void **) &(*object_dbl_ptr)->sequence_control_set_ptr,
        (void *)&scsInitData);
    if (return_error == EB_ErrorInsufficientResources) {
        return EB_ErrorInsufficientResources;
    }

    EB_CREATEMUTEX(EbHandle*, (*object_dbl_ptr)->config_mutex, sizeof(EbHandle), EB_MUTEX);


    return EB_ErrorNone;
}

extern EbErrorType sb_params_ctor(
    SequenceControlSet_t *sequence_control_set_ptr) {

    EbErrorType return_error = EB_ErrorNone;

    EB_MALLOC(SbParams_t*, sequence_control_set_ptr->sb_params_array, sizeof(SbParams_t) * ((MAX_PICTURE_WIDTH_SIZE + sequence_control_set_ptr->sb_sz - 1) / sequence_control_set_ptr->sb_sz) * ((MAX_PICTURE_HEIGHT_SIZE + sequence_control_set_ptr->sb_sz - 1) / sequence_control_set_ptr->sb_sz), EB_N_PTR);
    return return_error;
}

extern EbErrorType sb_params_init(
    SequenceControlSet_t *sequence_control_set_ptr) {

    EbErrorType return_error = EB_ErrorNone;
    uint16_t    sb_index;
    uint16_t    rasterScanCuIndex;
    uint16_t    md_scan_block_index;
    uint8_t   pictureLcuWidth = (uint8_t)((sequence_control_set_ptr->luma_width + sequence_control_set_ptr->sb_sz - 1) / sequence_control_set_ptr->sb_sz);
    uint8_t    pictureLcuHeight = (uint8_t)((sequence_control_set_ptr->luma_height + sequence_control_set_ptr->sb_sz - 1) / sequence_control_set_ptr->sb_sz);
    EB_MALLOC(SbParams_t*, sequence_control_set_ptr->sb_params_array, sizeof(SbParams_t) * pictureLcuWidth * pictureLcuHeight, EB_N_PTR);

    for (sb_index = 0; sb_index < pictureLcuWidth * pictureLcuHeight; ++sb_index) {
        sequence_control_set_ptr->sb_params_array[sb_index].horizontal_index = (uint8_t)(sb_index % pictureLcuWidth);
        sequence_control_set_ptr->sb_params_array[sb_index].vertical_index = (uint8_t)(sb_index / pictureLcuWidth);
        sequence_control_set_ptr->sb_params_array[sb_index].origin_x = sequence_control_set_ptr->sb_params_array[sb_index].horizontal_index * sequence_control_set_ptr->sb_sz;
        sequence_control_set_ptr->sb_params_array[sb_index].origin_y = sequence_control_set_ptr->sb_params_array[sb_index].vertical_index * sequence_control_set_ptr->sb_sz;

        sequence_control_set_ptr->sb_params_array[sb_index].width = (uint8_t)(((sequence_control_set_ptr->luma_width - sequence_control_set_ptr->sb_params_array[sb_index].origin_x) < sequence_control_set_ptr->sb_sz) ?
            sequence_control_set_ptr->luma_width - sequence_control_set_ptr->sb_params_array[sb_index].origin_x :
            sequence_control_set_ptr->sb_sz);

        sequence_control_set_ptr->sb_params_array[sb_index].height = (uint8_t)(((sequence_control_set_ptr->luma_height - sequence_control_set_ptr->sb_params_array[sb_index].origin_y) < sequence_control_set_ptr->sb_sz) ?
            sequence_control_set_ptr->luma_height - sequence_control_set_ptr->sb_params_array[sb_index].origin_y :
            sequence_control_set_ptr->sb_sz);

        sequence_control_set_ptr->sb_params_array[sb_index].is_complete_sb = (uint8_t)(((sequence_control_set_ptr->sb_params_array[sb_index].width == sequence_control_set_ptr->sb_sz) && (sequence_control_set_ptr->sb_params_array[sb_index].height == sequence_control_set_ptr->sb_sz)) ?
            1 :
            0);

        sequence_control_set_ptr->sb_params_array[sb_index].is_edge_sb = (sequence_control_set_ptr->sb_params_array[sb_index].origin_x < sequence_control_set_ptr->sb_sz) ||
            (sequence_control_set_ptr->sb_params_array[sb_index].origin_y < sequence_control_set_ptr->sb_sz) ||
            (sequence_control_set_ptr->sb_params_array[sb_index].origin_x > sequence_control_set_ptr->luma_width - sequence_control_set_ptr->sb_sz) ||
            (sequence_control_set_ptr->sb_params_array[sb_index].origin_y > sequence_control_set_ptr->luma_height - sequence_control_set_ptr->sb_sz) ? 1 : 0;


        uint8_t potential_logo_sb = 0;

        // 4K
        /*__ 14 __         __ 14__*/
        ///////////////////////////
        //         |          |         //
        //         8          8         //
        //___14_ |          |_14___//
        //                         //
        //                         //
        //-----------------------// |
        //                         // 8
        ///////////////////////////    |

        // 1080p/720P
        /*__ 7 __          __ 7__*/
        ///////////////////////////
        //         |          |         //
        //         4          4         //
        //___7__ |          |__7___//
        //                         //
        //                         //
        //-----------------------// |
        //                         // 4
        ///////////////////////////    |

        // 480P
        /*__ 3 __          __ 3__*/
        ///////////////////////////
        //         |          |         //
        //         2          2         //
        //___3__ |          |__3___//
        //                         //
        //                         //
        //-----------------------// |
        //                         // 2
        ///////////////////////////    |
        if (sequence_control_set_ptr->input_resolution <= INPUT_SIZE_576p_RANGE_OR_LOWER) {
            potential_logo_sb = (sequence_control_set_ptr->sb_params_array[sb_index].origin_x >= (sequence_control_set_ptr->luma_width - (3 * sequence_control_set_ptr->sb_sz))) && (sequence_control_set_ptr->sb_params_array[sb_index].origin_y < 2 * sequence_control_set_ptr->sb_sz) ? 1 : potential_logo_sb;
            potential_logo_sb = (sequence_control_set_ptr->sb_params_array[sb_index].origin_x < ((3 * sequence_control_set_ptr->sb_sz))) && (sequence_control_set_ptr->sb_params_array[sb_index].origin_y < 2 * sequence_control_set_ptr->sb_sz) ? 1 : potential_logo_sb;
            potential_logo_sb = (sequence_control_set_ptr->sb_params_array[sb_index].origin_y >= (sequence_control_set_ptr->luma_height - (2 * sequence_control_set_ptr->sb_sz))) ? 1 : potential_logo_sb;

        }
        else
            if (sequence_control_set_ptr->input_resolution < INPUT_SIZE_4K_RANGE) {
                potential_logo_sb = (sequence_control_set_ptr->sb_params_array[sb_index].origin_x >= (sequence_control_set_ptr->luma_width - (7 * sequence_control_set_ptr->sb_sz))) && (sequence_control_set_ptr->sb_params_array[sb_index].origin_y < 4 * sequence_control_set_ptr->sb_sz) ? 1 : potential_logo_sb;
                potential_logo_sb = (sequence_control_set_ptr->sb_params_array[sb_index].origin_x < ((7 * sequence_control_set_ptr->sb_sz))) && (sequence_control_set_ptr->sb_params_array[sb_index].origin_y < 4 * sequence_control_set_ptr->sb_sz) ? 1 : potential_logo_sb;
                potential_logo_sb = (sequence_control_set_ptr->sb_params_array[sb_index].origin_y >= (sequence_control_set_ptr->luma_height - (4 * sequence_control_set_ptr->sb_sz))) ? 1 : potential_logo_sb;

            }
            else {

                potential_logo_sb = (sequence_control_set_ptr->sb_params_array[sb_index].origin_x >= (sequence_control_set_ptr->luma_width - (14 * sequence_control_set_ptr->sb_sz))) && (sequence_control_set_ptr->sb_params_array[sb_index].origin_y < 8 * sequence_control_set_ptr->sb_sz) ? 1 : potential_logo_sb;
                potential_logo_sb = (sequence_control_set_ptr->sb_params_array[sb_index].origin_x < ((14 * sequence_control_set_ptr->sb_sz))) && (sequence_control_set_ptr->sb_params_array[sb_index].origin_y < 8 * sequence_control_set_ptr->sb_sz) ? 1 : potential_logo_sb;
                potential_logo_sb = (sequence_control_set_ptr->sb_params_array[sb_index].origin_y >= (sequence_control_set_ptr->luma_height - (8 * sequence_control_set_ptr->sb_sz))) ? 1 : potential_logo_sb;

            }
        sequence_control_set_ptr->sb_params_array[sb_index].potential_logo_sb = potential_logo_sb;

        for (rasterScanCuIndex = RASTER_SCAN_CU_INDEX_64x64; rasterScanCuIndex <= RASTER_SCAN_CU_INDEX_8x8_63; rasterScanCuIndex++) {

            sequence_control_set_ptr->sb_params_array[sb_index].raster_scan_cu_validity[rasterScanCuIndex] = ((sequence_control_set_ptr->sb_params_array[sb_index].origin_x + RASTER_SCAN_CU_X[rasterScanCuIndex] + RASTER_SCAN_CU_SIZE[rasterScanCuIndex] > sequence_control_set_ptr->luma_width) || (sequence_control_set_ptr->sb_params_array[sb_index].origin_y + RASTER_SCAN_CU_Y[rasterScanCuIndex] + RASTER_SCAN_CU_SIZE[rasterScanCuIndex] > sequence_control_set_ptr->luma_height)) ?
                EB_FALSE :
                EB_TRUE;
        }

         uint16_t max_block_count = sequence_control_set_ptr->max_block_cnt;

        for (md_scan_block_index = 0; md_scan_block_index < max_block_count; md_scan_block_index++) {

            const BlockGeom * blk_geom = Get_blk_geom_mds(md_scan_block_index);

            if (blk_geom->shape != PART_N)
                blk_geom = Get_blk_geom_mds(blk_geom->sqi_mds);

            sequence_control_set_ptr->sb_params_array[sb_index].block_is_inside_md_scan[md_scan_block_index] =
                ((sequence_control_set_ptr->sb_params_array[sb_index].origin_x + blk_geom->origin_x + blk_geom->bwidth > sequence_control_set_ptr->luma_width) ||
                (sequence_control_set_ptr->sb_params_array[sb_index].origin_y + blk_geom->origin_y + blk_geom->bheight > sequence_control_set_ptr->luma_height)) ?
                EB_FALSE :
                EB_TRUE;
        }

    }

    sequence_control_set_ptr->picture_width_in_sb = pictureLcuWidth;
    sequence_control_set_ptr->picture_height_in_sb = pictureLcuHeight;
    sequence_control_set_ptr->sb_total_count = pictureLcuWidth * pictureLcuHeight;

    return return_error;
}

EbErrorType sb_geom_init(SequenceControlSet_t * sequence_control_set_ptr)
{
    uint16_t    sb_index;
    uint16_t    md_scan_block_index;
    uint16_t   pictureLcuWidth = (sequence_control_set_ptr->luma_width + sequence_control_set_ptr->sb_size_pix - 1) / sequence_control_set_ptr->sb_size_pix;
    uint16_t    pictureLcuHeight = (sequence_control_set_ptr->luma_height + sequence_control_set_ptr->sb_size_pix - 1) / sequence_control_set_ptr->sb_size_pix;

    EB_MALLOC(SbGeom_t*, sequence_control_set_ptr->sb_geom, sizeof(SbGeom_t) * pictureLcuWidth * pictureLcuHeight, EB_N_PTR);


    for (sb_index = 0; sb_index < pictureLcuWidth * pictureLcuHeight; ++sb_index) {

        sequence_control_set_ptr->sb_geom[sb_index].horizontal_index = sb_index % pictureLcuWidth;
        sequence_control_set_ptr->sb_geom[sb_index].vertical_index = sb_index / pictureLcuWidth;
        sequence_control_set_ptr->sb_geom[sb_index].origin_x = sequence_control_set_ptr->sb_geom[sb_index].horizontal_index * sequence_control_set_ptr->sb_size_pix;
        sequence_control_set_ptr->sb_geom[sb_index].origin_y = sequence_control_set_ptr->sb_geom[sb_index].vertical_index * sequence_control_set_ptr->sb_size_pix;

        sequence_control_set_ptr->sb_geom[sb_index].width = (uint8_t)(((sequence_control_set_ptr->luma_width - sequence_control_set_ptr->sb_geom[sb_index].origin_x) < sequence_control_set_ptr->sb_size_pix) ?
            sequence_control_set_ptr->luma_width - sequence_control_set_ptr->sb_geom[sb_index].origin_x :
            sequence_control_set_ptr->sb_size_pix);

        sequence_control_set_ptr->sb_geom[sb_index].height = (uint8_t)(((sequence_control_set_ptr->luma_height - sequence_control_set_ptr->sb_geom[sb_index].origin_y) < sequence_control_set_ptr->sb_size_pix) ?
            sequence_control_set_ptr->luma_height - sequence_control_set_ptr->sb_geom[sb_index].origin_y :
            sequence_control_set_ptr->sb_size_pix);

        sequence_control_set_ptr->sb_geom[sb_index].is_complete_sb =
            (uint8_t)(((sequence_control_set_ptr->sb_geom[sb_index].width == sequence_control_set_ptr->sb_size_pix) && (sequence_control_set_ptr->sb_geom[sb_index].height == sequence_control_set_ptr->sb_size_pix)) ?
                1 :
                0);

        
        uint16_t max_block_count = sequence_control_set_ptr->max_block_cnt;

        for (md_scan_block_index = 0; md_scan_block_index < max_block_count ; md_scan_block_index++) {

            const BlockGeom * blk_geom = Get_blk_geom_mds(md_scan_block_index);

            if (blk_geom->shape != PART_N)
                blk_geom = Get_blk_geom_mds(blk_geom->sqi_mds);

            sequence_control_set_ptr->sb_geom[sb_index].block_is_inside_md_scan[md_scan_block_index] =
                ((sequence_control_set_ptr->sb_geom[sb_index].origin_x + blk_geom->origin_x + blk_geom->bwidth > sequence_control_set_ptr->luma_width) ||
                (sequence_control_set_ptr->sb_geom[sb_index].origin_y + blk_geom->origin_y + blk_geom->bheight > sequence_control_set_ptr->luma_height)) ?
                EB_FALSE :
                EB_TRUE;
        }

    }

    sequence_control_set_ptr->sb_tot_cnt = pictureLcuWidth * pictureLcuHeight;

    return 0;
}

