// clang-format off
/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>

#include "EbSequenceControlSet.h"
#include "EbUtility.h"

static void eb_sequence_control_set_dctor(EbPtr p)
{
    SequenceControlSet *obj = (SequenceControlSet*)p;
    EB_FREE_ARRAY(obj->sb_params_array);
    EB_FREE_ARRAY(obj->sb_geom);
}

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
    SequenceControlSet *sequence_control_set_ptr,
    EbPtr object_init_data_ptr)
{
    EbSequenceControlSetInitData *scsInitData = (EbSequenceControlSetInitData*)object_init_data_ptr;
    uint32_t segment_index;

    sequence_control_set_ptr->dctor = eb_sequence_control_set_dctor;

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
    if (scsInitData != EB_NULL)
        sequence_control_set_ptr->encode_context_ptr = scsInitData->encode_context_ptr;

    // Profile & ID
    sequence_control_set_ptr->chroma_format_idc = EB_YUV420;
    sequence_control_set_ptr->max_temporal_layers = 1;

    sequence_control_set_ptr->bits_for_picture_order_count = 16;
    sequence_control_set_ptr->subsampling_y = 1;
    sequence_control_set_ptr->subsampling_x = 1;

    sequence_control_set_ptr->encoder_bit_depth = 8;

    // Bitdepth
    //sequence_control_set_ptr->input_bitdepth = EB_8BIT;
    //sequence_control_set_ptr->output_bitdepth = EB_8BIT;

    // GOP Structure
    sequence_control_set_ptr->max_ref_count = 1;

    // LCU
    sequence_control_set_ptr->sb_sz = 64;
    sequence_control_set_ptr->max_sb_depth = 3;

    // CU
    sequence_control_set_ptr->max_cu_size = 64;
    sequence_control_set_ptr->min_cu_size = 8;
    sequence_control_set_ptr->max_intra_size = 32;
    sequence_control_set_ptr->min_intra_size = 8;
    sequence_control_set_ptr->intra4x4_flag = EB_TRUE;
    // Rate Control
    sequence_control_set_ptr->rate_control_mode = 0;
    sequence_control_set_ptr->target_bitrate = 0x1000;
    sequence_control_set_ptr->available_bandwidth = 0x1000;

    // Quantization
    sequence_control_set_ptr->static_config.qp = 20;
    // Initialize SB params
    EB_MALLOC_ARRAY(sequence_control_set_ptr->sb_params_array,
        ((MAX_PICTURE_WIDTH_SIZE + sequence_control_set_ptr->sb_sz - 1) / sequence_control_set_ptr->sb_sz) *
        ((MAX_PICTURE_HEIGHT_SIZE + sequence_control_set_ptr->sb_sz - 1) / sequence_control_set_ptr->sb_sz));

    sequence_control_set_ptr->seq_header.frame_width_bits = 16;
    sequence_control_set_ptr->seq_header.frame_height_bits = 16;
    sequence_control_set_ptr->seq_header.frame_id_numbers_present_flag = 0;
    //    cm->large_scale_tile ? 0 : cm->error_resilient_mode;
    sequence_control_set_ptr->seq_header.frame_id_length = FRAME_ID_LENGTH;
    sequence_control_set_ptr->seq_header.delta_frame_id_length = DELTA_FRAME_ID_LENGTH;

    if (scsInitData && scsInitData->sb_size == 128)
    {
        sequence_control_set_ptr->seq_header.sb_size = BLOCK_128X128;
        sequence_control_set_ptr->sb_size_pix = 128;
        sequence_control_set_ptr->max_block_cnt = 4421;

        sequence_control_set_ptr->seq_header.sb_mi_size = 32;        // Size of the superblock in units of MI blocks
        sequence_control_set_ptr->seq_header.sb_size_log2 = 5;
    }
    else
    {
        sequence_control_set_ptr->seq_header.sb_size = BLOCK_64X64;
        sequence_control_set_ptr->sb_size_pix = 64;
        sequence_control_set_ptr->max_block_cnt = 1101;

        sequence_control_set_ptr->seq_header.sb_mi_size = 16;        // Size of the superblock in units of MI blocks
        sequence_control_set_ptr->seq_header.sb_size_log2 = 4;
    }
    // 0 - disable dual interpolation filter
    // 1 - enable vertical and horiz filter selection
    sequence_control_set_ptr->seq_header.enable_dual_filter = 0;
    sequence_control_set_ptr->seq_header.order_hint_info.enable_order_hint = 1;
    // 0 - disable order hint, and related tools:
    // jnt_comp, ref_frame_mvs, frame_sign_bias
    // if 0, enable_jnt_comp must be set zs 0.
    sequence_control_set_ptr->seq_header.order_hint_info.enable_jnt_comp = 0;

    sequence_control_set_ptr->seq_header.order_hint_info.order_hint_bits = sequence_control_set_ptr->seq_header.order_hint_info.enable_order_hint ? (6+1) : (-1+1);

    sequence_control_set_ptr->seq_header.seq_force_screen_content_tools = 2;
    // 0 - force off
    // 1 - force on
    // 2 - adaptive
    sequence_control_set_ptr->seq_header.seq_force_integer_mv = 2;  // 0 - Not to force. MV can be in 1/4 or 1/8
    // 1 - force to integer
    // 2 - adaptive

    sequence_control_set_ptr->seq_header.order_hint_info.enable_ref_frame_mvs = 1;
#if NO_ENCDEC || SHUT_FILTERING
    sequence_control_set_ptr->seq_header.enable_cdef = 0;

    if (sequence_control_set_ptr->static_config.enable_restoration_filtering == DEFAULT)
        sequence_control_set_ptr->seq_header.enable_restoration = 0;
    else
        sequence_control_set_ptr->seq_header.enable_restoration = (uint8_t)sequence_control_set_ptr->static_config.enable_restoration_filtering;
#else
    sequence_control_set_ptr->seq_header.enable_cdef = 1;
    if (sequence_control_set_ptr->static_config.enable_restoration_filtering == DEFAULT)
        sequence_control_set_ptr->seq_header.enable_restoration = 1;
    else
        sequence_control_set_ptr->seq_header.enable_restoration = (uint8_t)sequence_control_set_ptr->static_config.enable_restoration_filtering;
#endif

    sequence_control_set_ptr->film_grain_random_seed = 7391;
    sequence_control_set_ptr->reference_count = 4;

    return EB_ErrorNone;
}

EbErrorType eb_sequence_control_set_creator(
    EbPtr *object_dbl_ptr,
    EbPtr object_init_data_ptr)
{
    SequenceControlSet* obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, eb_sequence_control_set_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;
    return EB_ErrorNone;
}

/************************************************
 * Sequence Control Set Copy
 ************************************************/
EbErrorType copy_sequence_control_set(
    SequenceControlSet *dst,
    SequenceControlSet *src)
{
    uint32_t  writeCount = 0;

    dst->static_config = src->static_config;                            writeCount += sizeof(EbSvtAv1EncConfiguration);
    dst->encode_context_ptr = src->encode_context_ptr;                        writeCount += sizeof(EncodeContext*);
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
    dst->seq_header.max_frame_width = src->seq_header.max_frame_width;                               writeCount += sizeof(uint32_t);
    dst->seq_header.max_frame_height = src->seq_header.max_frame_height;                              writeCount += sizeof(uint32_t);
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
    //dst->input_bitdepth = src->input_bitdepth;                           writeCount += sizeof(EB_BITDEPTH);
    //dst->output_bitdepth = src->output_bitdepth;                          writeCount += sizeof(EB_BITDEPTH);
    dst->encoder_bit_depth = src->encoder_bit_depth;                      writeCount += sizeof(uint32_t);
    dst->subsampling_x = src->subsampling_x;                writeCount += sizeof(uint16_t);
    dst->subsampling_y = src->subsampling_y;                writeCount += sizeof(uint16_t);
    dst->pred_struct_ptr = src->pred_struct_ptr;                           writeCount += sizeof(PredictionStructure*);
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
    dst->target_bitrate = src->target_bitrate;                           writeCount += sizeof(uint32_t);
    dst->available_bandwidth = src->available_bandwidth;                      writeCount += sizeof(uint32_t);
    dst->static_config.qp = src->static_config.qp;                                      writeCount += sizeof(uint32_t);
    dst->film_grain_denoise_strength = src->film_grain_denoise_strength;          writeCount += sizeof(int32_t);
    dst->seq_header.film_grain_params_present = src->seq_header.film_grain_params_present;              writeCount += sizeof(int32_t);
    dst->seq_header.film_grain_params_present = src->seq_header.film_grain_params_present;              writeCount += sizeof(int32_t);
    dst->picture_control_set_pool_init_count = src->picture_control_set_pool_init_count;            writeCount += sizeof(int32_t);
    dst->picture_control_set_pool_init_count_child = src->picture_control_set_pool_init_count_child; writeCount += sizeof(int32_t);
    dst->pa_reference_picture_buffer_init_count = src->pa_reference_picture_buffer_init_count; writeCount += sizeof(int32_t);
    dst->reference_picture_buffer_init_count = src->reference_picture_buffer_init_count; writeCount += sizeof(int32_t);
    dst->input_buffer_fifo_init_count = src->input_buffer_fifo_init_count; writeCount += sizeof(int32_t);
    dst->overlay_input_picture_buffer_init_count = src->overlay_input_picture_buffer_init_count; writeCount += sizeof(int32_t);

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
    dst->reference_count = src->reference_count; writeCount += sizeof(uint32_t);
    for (uint8_t i = 0; i< MAX_HIERARCHICAL_LEVEL; i++) {
        dst->me_segment_column_count_array[i] = src->me_segment_column_count_array[i];
        dst->me_segment_row_count_array[i] = src->me_segment_row_count_array[i];
        dst->enc_dec_segment_col_count_array[i] = src->enc_dec_segment_col_count_array[i];
        dst->enc_dec_segment_row_count_array[i] = src->enc_dec_segment_row_count_array[i];
    }

    dst->cdef_segment_column_count = src->cdef_segment_column_count;
    dst->cdef_segment_row_count = src->cdef_segment_row_count;

    dst->rest_segment_column_count = src->rest_segment_column_count;
    dst->rest_segment_row_count = src->rest_segment_row_count;
    dst->mrp_mode       = src->mrp_mode;
    dst->nsq_present    = src->nsq_present;
    dst->cdf_mode       = src->cdf_mode;
    dst->down_sampling_method_me_search = src->down_sampling_method_me_search;
    dst->tf_segment_column_count = src->tf_segment_column_count;
    dst->tf_segment_row_count = src->tf_segment_row_count;
    dst->over_boundary_block_mode = src->over_boundary_block_mode;
    dst->mfmv_enabled = src->mfmv_enabled;
#if TWO_PASS
    dst->use_input_stat_file = src->use_input_stat_file;
    dst->use_output_stat_file = src->use_output_stat_file;
#endif
    dst->scd_delay = src->scd_delay;
    return EB_ErrorNone;
}

extern EbErrorType derive_input_resolution(
    SequenceControlSet *sequenceControlSetPtr,
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

static void eb_sequence_control_set_instance_dctor(EbPtr p)
{
   EbSequenceControlSetInstance* obj = (EbSequenceControlSetInstance*)p;
   EB_DELETE(obj->encode_context_ptr);
   EB_DELETE(obj->sequence_control_set_ptr);
   EB_DESTROY_MUTEX(obj->config_mutex);
}

EbErrorType eb_sequence_control_set_instance_ctor(
    EbSequenceControlSetInstance *object_dbl_ptr)
{
    EbSequenceControlSetInitData scsInitData;

    object_dbl_ptr->dctor = eb_sequence_control_set_instance_dctor;

    EB_NEW(object_dbl_ptr->encode_context_ptr, encode_context_ctor, EB_NULL);
    scsInitData.encode_context_ptr = object_dbl_ptr->encode_context_ptr;

    scsInitData.sb_size = 64;

    EB_NEW(object_dbl_ptr->sequence_control_set_ptr, eb_sequence_control_set_ctor, (void *)&scsInitData);
    EB_CREATE_MUTEX(object_dbl_ptr->config_mutex);

    return EB_ErrorNone;
}

extern EbErrorType sb_params_init(
    SequenceControlSet *sequence_control_set_ptr) {
    EbErrorType return_error = EB_ErrorNone;
    uint16_t    sb_index;
    uint16_t    rasterScanCuIndex;
    uint8_t   pictureLcuWidth = (uint8_t)((sequence_control_set_ptr->seq_header.max_frame_width + sequence_control_set_ptr->sb_sz - 1) / sequence_control_set_ptr->sb_sz);
    uint8_t    pictureLcuHeight = (uint8_t)((sequence_control_set_ptr->seq_header.max_frame_height + sequence_control_set_ptr->sb_sz - 1) / sequence_control_set_ptr->sb_sz);
    //free old one;
    EB_FREE_ARRAY(sequence_control_set_ptr->sb_params_array);

    EB_MALLOC_ARRAY(sequence_control_set_ptr->sb_params_array, pictureLcuWidth * pictureLcuHeight);

    for (sb_index = 0; sb_index < pictureLcuWidth * pictureLcuHeight; ++sb_index) {
        sequence_control_set_ptr->sb_params_array[sb_index].horizontal_index = (uint8_t)(sb_index % pictureLcuWidth);
        sequence_control_set_ptr->sb_params_array[sb_index].vertical_index = (uint8_t)(sb_index / pictureLcuWidth);
        sequence_control_set_ptr->sb_params_array[sb_index].origin_x = sequence_control_set_ptr->sb_params_array[sb_index].horizontal_index * sequence_control_set_ptr->sb_sz;
        sequence_control_set_ptr->sb_params_array[sb_index].origin_y = sequence_control_set_ptr->sb_params_array[sb_index].vertical_index * sequence_control_set_ptr->sb_sz;

        sequence_control_set_ptr->sb_params_array[sb_index].width = (uint8_t)(((sequence_control_set_ptr->seq_header.max_frame_width - sequence_control_set_ptr->sb_params_array[sb_index].origin_x) < sequence_control_set_ptr->sb_sz) ?
            sequence_control_set_ptr->seq_header.max_frame_width - sequence_control_set_ptr->sb_params_array[sb_index].origin_x :
            sequence_control_set_ptr->sb_sz);

        sequence_control_set_ptr->sb_params_array[sb_index].height = (uint8_t)(((sequence_control_set_ptr->seq_header.max_frame_height - sequence_control_set_ptr->sb_params_array[sb_index].origin_y) < sequence_control_set_ptr->sb_sz) ?
            sequence_control_set_ptr->seq_header.max_frame_height - sequence_control_set_ptr->sb_params_array[sb_index].origin_y :
            sequence_control_set_ptr->sb_sz);

        sequence_control_set_ptr->sb_params_array[sb_index].is_complete_sb = (uint8_t)(((sequence_control_set_ptr->sb_params_array[sb_index].width == sequence_control_set_ptr->sb_sz) && (sequence_control_set_ptr->sb_params_array[sb_index].height == sequence_control_set_ptr->sb_sz)) ?
            1 :
            0);

        sequence_control_set_ptr->sb_params_array[sb_index].is_edge_sb = (sequence_control_set_ptr->sb_params_array[sb_index].origin_x < sequence_control_set_ptr->sb_sz) ||
            (sequence_control_set_ptr->sb_params_array[sb_index].origin_y < sequence_control_set_ptr->sb_sz) ||
            (sequence_control_set_ptr->sb_params_array[sb_index].origin_x > sequence_control_set_ptr->seq_header.max_frame_width - sequence_control_set_ptr->sb_sz) ||
            (sequence_control_set_ptr->sb_params_array[sb_index].origin_y > sequence_control_set_ptr->seq_header.max_frame_height - sequence_control_set_ptr->sb_sz) ? 1 : 0;

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
            potential_logo_sb = (sequence_control_set_ptr->sb_params_array[sb_index].origin_x >= (sequence_control_set_ptr->seq_header.max_frame_width - (3 * sequence_control_set_ptr->sb_sz))) && (sequence_control_set_ptr->sb_params_array[sb_index].origin_y < 2 * sequence_control_set_ptr->sb_sz) ? 1 : potential_logo_sb;
            potential_logo_sb = (sequence_control_set_ptr->sb_params_array[sb_index].origin_x < ((3 * sequence_control_set_ptr->sb_sz))) && (sequence_control_set_ptr->sb_params_array[sb_index].origin_y < 2 * sequence_control_set_ptr->sb_sz) ? 1 : potential_logo_sb;
            potential_logo_sb = (sequence_control_set_ptr->sb_params_array[sb_index].origin_y >= (sequence_control_set_ptr->seq_header.max_frame_height - (2 * sequence_control_set_ptr->sb_sz))) ? 1 : potential_logo_sb;
        }
        else
            if (sequence_control_set_ptr->input_resolution < INPUT_SIZE_4K_RANGE) {
                potential_logo_sb = (sequence_control_set_ptr->sb_params_array[sb_index].origin_x >= (sequence_control_set_ptr->seq_header.max_frame_width - (7 * sequence_control_set_ptr->sb_sz))) && (sequence_control_set_ptr->sb_params_array[sb_index].origin_y < 4 * sequence_control_set_ptr->sb_sz) ? 1 : potential_logo_sb;
                potential_logo_sb = (sequence_control_set_ptr->sb_params_array[sb_index].origin_x < ((7 * sequence_control_set_ptr->sb_sz))) && (sequence_control_set_ptr->sb_params_array[sb_index].origin_y < 4 * sequence_control_set_ptr->sb_sz) ? 1 : potential_logo_sb;
                potential_logo_sb = (sequence_control_set_ptr->sb_params_array[sb_index].origin_y >= (sequence_control_set_ptr->seq_header.max_frame_height - (4 * sequence_control_set_ptr->sb_sz))) ? 1 : potential_logo_sb;
            }
            else {
                potential_logo_sb = (sequence_control_set_ptr->sb_params_array[sb_index].origin_x >= (sequence_control_set_ptr->seq_header.max_frame_width - (14 * sequence_control_set_ptr->sb_sz))) && (sequence_control_set_ptr->sb_params_array[sb_index].origin_y < 8 * sequence_control_set_ptr->sb_sz) ? 1 : potential_logo_sb;
                potential_logo_sb = (sequence_control_set_ptr->sb_params_array[sb_index].origin_x < ((14 * sequence_control_set_ptr->sb_sz))) && (sequence_control_set_ptr->sb_params_array[sb_index].origin_y < 8 * sequence_control_set_ptr->sb_sz) ? 1 : potential_logo_sb;
                potential_logo_sb = (sequence_control_set_ptr->sb_params_array[sb_index].origin_y >= (sequence_control_set_ptr->seq_header.max_frame_height - (8 * sequence_control_set_ptr->sb_sz))) ? 1 : potential_logo_sb;
            }
        sequence_control_set_ptr->sb_params_array[sb_index].potential_logo_sb = potential_logo_sb;

        for (rasterScanCuIndex = RASTER_SCAN_CU_INDEX_64x64; rasterScanCuIndex <= RASTER_SCAN_CU_INDEX_8x8_63; rasterScanCuIndex++) {
            sequence_control_set_ptr->sb_params_array[sb_index].raster_scan_cu_validity[rasterScanCuIndex] = ((sequence_control_set_ptr->sb_params_array[sb_index].origin_x + raster_scan_cu_x[rasterScanCuIndex] + raster_scan_cu_size[rasterScanCuIndex] > sequence_control_set_ptr->seq_header.max_frame_width) || (sequence_control_set_ptr->sb_params_array[sb_index].origin_y + raster_scan_cu_y[rasterScanCuIndex] + raster_scan_cu_size[rasterScanCuIndex] > sequence_control_set_ptr->seq_header.max_frame_height)) ?
                EB_FALSE :
                EB_TRUE;
        }
    }

    sequence_control_set_ptr->picture_width_in_sb = pictureLcuWidth;
    sequence_control_set_ptr->picture_height_in_sb = pictureLcuHeight;
    sequence_control_set_ptr->sb_total_count = pictureLcuWidth * pictureLcuHeight;

    return return_error;
}

EbErrorType sb_geom_init(SequenceControlSet * sequence_control_set_ptr)
{
    uint16_t    sb_index;
    uint16_t    md_scan_block_index;
    uint16_t   pictureLcuWidth = (sequence_control_set_ptr->seq_header.max_frame_width + sequence_control_set_ptr->sb_size_pix - 1) / sequence_control_set_ptr->sb_size_pix;
    uint16_t    pictureLcuHeight = (sequence_control_set_ptr->seq_header.max_frame_height + sequence_control_set_ptr->sb_size_pix - 1) / sequence_control_set_ptr->sb_size_pix;


    EB_FREE_ARRAY(sequence_control_set_ptr->sb_geom);
    EB_MALLOC_ARRAY(sequence_control_set_ptr->sb_geom, pictureLcuWidth * pictureLcuHeight);

    for (sb_index = 0; sb_index < pictureLcuWidth * pictureLcuHeight; ++sb_index) {
        sequence_control_set_ptr->sb_geom[sb_index].horizontal_index = sb_index % pictureLcuWidth;
        sequence_control_set_ptr->sb_geom[sb_index].vertical_index = sb_index / pictureLcuWidth;
        sequence_control_set_ptr->sb_geom[sb_index].origin_x = sequence_control_set_ptr->sb_geom[sb_index].horizontal_index * sequence_control_set_ptr->sb_size_pix;
        sequence_control_set_ptr->sb_geom[sb_index].origin_y = sequence_control_set_ptr->sb_geom[sb_index].vertical_index * sequence_control_set_ptr->sb_size_pix;

        sequence_control_set_ptr->sb_geom[sb_index].width = (uint8_t)(((sequence_control_set_ptr->seq_header.max_frame_width - sequence_control_set_ptr->sb_geom[sb_index].origin_x) < sequence_control_set_ptr->sb_size_pix) ?
            sequence_control_set_ptr->seq_header.max_frame_width - sequence_control_set_ptr->sb_geom[sb_index].origin_x :
            sequence_control_set_ptr->sb_size_pix);

        sequence_control_set_ptr->sb_geom[sb_index].height = (uint8_t)(((sequence_control_set_ptr->seq_header.max_frame_height - sequence_control_set_ptr->sb_geom[sb_index].origin_y) < sequence_control_set_ptr->sb_size_pix) ?
            sequence_control_set_ptr->seq_header.max_frame_height - sequence_control_set_ptr->sb_geom[sb_index].origin_y :
            sequence_control_set_ptr->sb_size_pix);

        sequence_control_set_ptr->sb_geom[sb_index].is_complete_sb =
            (uint8_t)(((sequence_control_set_ptr->sb_geom[sb_index].width == sequence_control_set_ptr->sb_size_pix) && (sequence_control_set_ptr->sb_geom[sb_index].height == sequence_control_set_ptr->sb_size_pix)) ?
                1 :
                0);

        uint16_t max_block_count = sequence_control_set_ptr->max_block_cnt;

        for (md_scan_block_index = 0; md_scan_block_index < max_block_count ; md_scan_block_index++) {
            const BlockGeom * blk_geom = get_blk_geom_mds(md_scan_block_index);
            if (sequence_control_set_ptr->over_boundary_block_mode == 1) {
                sequence_control_set_ptr->sb_geom[sb_index].block_is_allowed[md_scan_block_index] =
                    ((sequence_control_set_ptr->sb_geom[sb_index].origin_x + blk_geom->origin_x + blk_geom->bwidth / 2 < sequence_control_set_ptr->seq_header.max_frame_width) &&
                    (sequence_control_set_ptr->sb_geom[sb_index].origin_y + blk_geom->origin_y + blk_geom->bheight / 2 < sequence_control_set_ptr->seq_header.max_frame_height)) ?
                    EB_TRUE :
                    EB_FALSE;

                // Temporary if the cropped width is not 4, 8, 16, 32, 64 and 128, the block is not allowed. To be removed after intrinsic functions for NxM spatial_full_distortion_kernel_func_ptr_array are added
                int32_t cropped_width = MIN(blk_geom->bwidth, sequence_control_set_ptr->seq_header.max_frame_width - (sequence_control_set_ptr->sb_geom[sb_index].origin_x + blk_geom->origin_x));
                if (cropped_width != 4 && cropped_width != 8 && cropped_width != 16 && cropped_width != 32 && cropped_width != 64 && cropped_width != 128)
                    sequence_control_set_ptr->sb_geom[sb_index].block_is_allowed[md_scan_block_index] = EB_FALSE;

                if (blk_geom->shape != PART_N)
                    blk_geom = get_blk_geom_mds(blk_geom->sqi_mds);
                sequence_control_set_ptr->sb_geom[sb_index].block_is_inside_md_scan[md_scan_block_index] =
                    ((sequence_control_set_ptr->sb_geom[sb_index].origin_x >= sequence_control_set_ptr->seq_header.max_frame_width) ||
                    (sequence_control_set_ptr->sb_geom[sb_index].origin_y >= sequence_control_set_ptr->seq_header.max_frame_height)) ?
                    EB_FALSE :
                    EB_TRUE;
            }
            else {
                if (blk_geom->shape != PART_N)
                    blk_geom = get_blk_geom_mds(blk_geom->sqi_mds);

                sequence_control_set_ptr->sb_geom[sb_index].block_is_allowed[md_scan_block_index] =
                    ((sequence_control_set_ptr->sb_geom[sb_index].origin_x + blk_geom->origin_x + blk_geom->bwidth > sequence_control_set_ptr->seq_header.max_frame_width) ||
                    (sequence_control_set_ptr->sb_geom[sb_index].origin_y + blk_geom->origin_y + blk_geom->bheight > sequence_control_set_ptr->seq_header.max_frame_height)) ?
                    EB_FALSE :
                    EB_TRUE;

                sequence_control_set_ptr->sb_geom[sb_index].block_is_inside_md_scan[md_scan_block_index] =
                    ((sequence_control_set_ptr->sb_geom[sb_index].origin_x + blk_geom->origin_x + blk_geom->bwidth > sequence_control_set_ptr->seq_header.max_frame_width) ||
                    (sequence_control_set_ptr->sb_geom[sb_index].origin_y + blk_geom->origin_y + blk_geom->bheight > sequence_control_set_ptr->seq_header.max_frame_height)) ?
                    EB_FALSE :
                    EB_TRUE;

            }
        }
    }

    sequence_control_set_ptr->sb_tot_cnt = pictureLcuWidth * pictureLcuHeight;

    return 0;
}
// clang-format on
