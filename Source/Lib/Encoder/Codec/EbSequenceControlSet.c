/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>

#include "EbSequenceControlSet.h"
#include "EbUtility.h"

static void eb_sequence_control_set_dctor(EbPtr p) {
    SequenceControlSet *obj = (SequenceControlSet *)p;
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
EbErrorType eb_sequence_control_set_ctor(SequenceControlSet *scs_ptr, EbPtr object_init_data_ptr) {
    EbSequenceControlSetInitData *scs_init_data =
        (EbSequenceControlSetInitData *)object_init_data_ptr;
    uint32_t segment_index;

    scs_ptr->dctor = eb_sequence_control_set_dctor;

    scs_ptr->static_config.sb_sz           = 64;
    scs_ptr->static_config.partition_depth = 4;
    scs_ptr->static_config.qp              = 32;

    // Segments
    for (segment_index = 0; segment_index < MAX_TEMPORAL_LAYERS; ++segment_index) {
        scs_ptr->me_segment_column_count_array[segment_index]   = 1;
        scs_ptr->me_segment_row_count_array[segment_index]      = 1;
        scs_ptr->enc_dec_segment_col_count_array[segment_index] = 1;
        scs_ptr->enc_dec_segment_row_count_array[segment_index] = 1;
    }

    // Encode Context
    if (scs_init_data != NULL) scs_ptr->encode_context_ptr = scs_init_data->encode_context_ptr;

    // Profile & ID
    scs_ptr->chroma_format_idc   = EB_YUV420;
    scs_ptr->max_temporal_layers = 1;

    scs_ptr->bits_for_picture_order_count = 16;
    scs_ptr->subsampling_y                = 1;
    scs_ptr->subsampling_x                = 1;

    scs_ptr->encoder_bit_depth = 8;

    // Bitdepth
    //scs_ptr->input_bitdepth = EB_8BIT;
    //scs_ptr->output_bitdepth = EB_8BIT;

    // GOP Structure
    scs_ptr->max_ref_count = 1;

    // SB
    scs_ptr->sb_sz        = 64;
    scs_ptr->max_sb_depth = 3;

    // CU
    scs_ptr->max_blk_size   = 64;
    scs_ptr->min_blk_size   = 8;
    scs_ptr->max_intra_size = 32;
    scs_ptr->min_intra_size = 8;
    // Rate Control
    scs_ptr->target_bitrate      = 0x1000;
    // Quantization
    scs_ptr->static_config.qp = 20;
    // Initialize SB params
    EB_MALLOC_ARRAY(scs_ptr->sb_params_array,
                    ((MAX_PICTURE_WIDTH_SIZE + scs_ptr->sb_sz - 1) / scs_ptr->sb_sz) *
                        ((MAX_PICTURE_HEIGHT_SIZE + scs_ptr->sb_sz - 1) / scs_ptr->sb_sz));

    scs_ptr->seq_header.frame_width_bits              = 16;
    scs_ptr->seq_header.frame_height_bits             = 16;
    scs_ptr->seq_header.frame_id_numbers_present_flag = 0;
    //    cm->large_scale_tile ? 0 : cm->error_resilient_mode;
    scs_ptr->seq_header.frame_id_length       = FRAME_ID_LENGTH;
    scs_ptr->seq_header.delta_frame_id_length = DELTA_FRAME_ID_LENGTH;

    if (scs_init_data && scs_init_data->sb_size == 128) {
        scs_ptr->seq_header.sb_size = BLOCK_128X128;
        scs_ptr->sb_size_pix        = 128;
        scs_ptr->max_block_cnt      = 4421;

        scs_ptr->seq_header.sb_mi_size   = 32; // Size of the superblock in units of MI blocks
        scs_ptr->seq_header.sb_size_log2 = 5;
    } else {
        scs_ptr->seq_header.sb_size = BLOCK_64X64;
        scs_ptr->sb_size_pix        = 64;
        scs_ptr->max_block_cnt      = 1101;

        scs_ptr->seq_header.sb_mi_size   = 16; // Size of the superblock in units of MI blocks
        scs_ptr->seq_header.sb_size_log2 = 4;
    }
    // 0 - disable dual interpolation filter
    // 1 - enable vertical and horiz filter selection
    scs_ptr->seq_header.enable_dual_filter                = 0;
    scs_ptr->seq_header.order_hint_info.enable_order_hint = 1;
    // 0 - disable order hint, and related tools:
    // jnt_comp, ref_frame_mvs, frame_sign_bias
    // if 0, enable_jnt_comp must be set zs 0.
    scs_ptr->seq_header.order_hint_info.enable_jnt_comp = 0;

    scs_ptr->seq_header.order_hint_info.order_hint_bits =
        scs_ptr->seq_header.order_hint_info.enable_order_hint ? (6 + 1) : (-1 + 1);

    scs_ptr->seq_header.seq_force_screen_content_tools = 2;
    // 0 - force off
    // 1 - force on
    // 2 - adaptive
    scs_ptr->seq_header.seq_force_integer_mv = 2; // 0 - Not to force. MV can be in 1/4 or 1/8
    // 1 - force to integer
    // 2 - adaptive

    scs_ptr->seq_header.order_hint_info.enable_ref_frame_mvs = 1;
#if NO_ENCDEC
    if (scs_ptr->static_config.cdef_mode == DEFAULT)
        scs_ptr->seq_header.enable_cdef = 0;
    else
        scs_ptr->seq_header.enable_cdef = (uint8_t)(scs_ptr->static_config.cdef_mode > 0);

    if (scs_ptr->static_config.enable_restoration_filtering == DEFAULT)
        scs_ptr->seq_header.enable_restoration = 0;
    else
        scs_ptr->seq_header.enable_restoration =
            (uint8_t)scs_ptr->static_config.enable_restoration_filtering;
#else
    if (scs_ptr->static_config.cdef_mode == DEFAULT)
        scs_ptr->seq_header.enable_cdef = 1;
    else
        scs_ptr->seq_header.enable_cdef = (uint8_t)(scs_ptr->static_config.cdef_mode > 0);

    if (scs_ptr->static_config.enable_restoration_filtering == DEFAULT)
        scs_ptr->seq_header.enable_restoration = 1;
    else
        scs_ptr->seq_header.enable_restoration =
            (uint8_t)scs_ptr->static_config.enable_restoration_filtering;
#endif

    if (scs_ptr->static_config.enable_intra_edge_filter == DEFAULT)
        scs_ptr->seq_header.enable_intra_edge_filter = 1;
    else
        scs_ptr->seq_header.enable_intra_edge_filter = (uint8_t)scs_ptr->static_config.enable_intra_edge_filter;

    if (scs_ptr->static_config.pic_based_rate_est == DEFAULT)
#if PIC_BASED_RE_OFF
        scs_ptr->seq_header.pic_based_rate_est = 0;
#else
        scs_ptr->seq_header.pic_based_rate_est = 1;
#endif
    else
        scs_ptr->seq_header.pic_based_rate_est = (uint8_t)scs_ptr->static_config.pic_based_rate_est;

    if (scs_ptr->static_config.enable_warped_motion == DEFAULT)
        scs_ptr->seq_header.enable_warped_motion = 1;
    else
        scs_ptr->seq_header.enable_warped_motion = (uint8_t)scs_ptr->static_config.enable_warped_motion;

    scs_ptr->film_grain_random_seed = 7391;
    scs_ptr->reference_count        = 4;

    return EB_ErrorNone;
}

EbErrorType eb_sequence_control_set_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr) {
    SequenceControlSet *obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, eb_sequence_control_set_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;
    return EB_ErrorNone;
}

/************************************************
 * Sequence Control Set Copy
 ************************************************/
EbErrorType copy_sequence_control_set(SequenceControlSet *dst, SequenceControlSet *src) {
    uint32_t write_count = 0;

    dst->static_config = src->static_config;
    write_count += sizeof(EbSvtAv1EncConfiguration);
    dst->encode_context_ptr = src->encode_context_ptr;
    write_count += sizeof(EncodeContext *);
    dst->chroma_format_idc = src->chroma_format_idc;
    write_count += sizeof(uint32_t);
    dst->max_temporal_layers = src->max_temporal_layers;
    write_count += sizeof(uint32_t);
    dst->bits_for_picture_order_count = src->bits_for_picture_order_count;
    write_count += sizeof(uint32_t);
    dst->max_input_luma_width = src->max_input_luma_width;
    write_count += sizeof(uint32_t);
    dst->max_input_luma_height = src->max_input_luma_height;
    write_count += sizeof(uint32_t);
    dst->max_input_chroma_height = src->max_input_chroma_height;
    write_count += sizeof(uint32_t);
    dst->max_input_chroma_width = src->max_input_chroma_width;
    write_count += sizeof(uint32_t);
    dst->max_input_pad_right = src->max_input_pad_right;
    write_count += sizeof(uint32_t);
    dst->max_input_pad_bottom = src->max_input_pad_bottom;
    write_count += sizeof(uint32_t);
    dst->seq_header.max_frame_width = src->seq_header.max_frame_width;
    write_count += sizeof(uint32_t);
    dst->seq_header.max_frame_height = src->seq_header.max_frame_height;
    write_count += sizeof(uint32_t);
    dst->chroma_width = src->chroma_width;
    write_count += sizeof(uint32_t);
    dst->chroma_height = src->chroma_height;
    write_count += sizeof(uint32_t);
    dst->pad_right = src->pad_right;
    write_count += sizeof(uint32_t);
    dst->pad_bottom = src->pad_bottom;
    write_count += sizeof(uint32_t);
    dst->frame_rate = src->frame_rate;
    write_count += sizeof(uint32_t);
    //dst->input_bitdepth = src->input_bitdepth;                           write_count += sizeof(EB_BITDEPTH);
    //dst->output_bitdepth = src->output_bitdepth;                          write_count += sizeof(EB_BITDEPTH);
    dst->encoder_bit_depth = src->encoder_bit_depth;
    write_count += sizeof(uint32_t);
    dst->subsampling_x = src->subsampling_x;
    write_count += sizeof(uint16_t);
    dst->subsampling_y = src->subsampling_y;
    write_count += sizeof(uint16_t);
    dst->pred_struct_ptr = src->pred_struct_ptr;
    write_count += sizeof(PredictionStructure *);
    dst->intra_period_length = src->intra_period_length;
    write_count += sizeof(int32_t);
    dst->intra_refresh_type = src->intra_refresh_type;
    write_count += sizeof(uint32_t);
    dst->max_ref_count = src->max_ref_count;
    write_count += sizeof(uint32_t);
    dst->sb_sz = src->sb_sz;
    write_count += sizeof(uint32_t);
    dst->max_sb_depth = src->max_sb_depth;
    write_count += sizeof(uint32_t);
    dst->max_blk_size = src->max_blk_size;
    write_count += sizeof(uint32_t);
    dst->min_blk_size = src->min_blk_size;
    write_count += sizeof(uint32_t);
    dst->max_intra_size = src->max_intra_size;
    write_count += sizeof(uint32_t);
    dst->min_intra_size = src->min_intra_size;
    write_count += sizeof(uint32_t);
    dst->target_bitrate = src->target_bitrate;
    write_count += sizeof(uint32_t);
    dst->static_config.qp = src->static_config.qp;
    write_count += sizeof(uint32_t);
    dst->film_grain_denoise_strength = src->film_grain_denoise_strength;
    write_count += sizeof(int32_t);
    dst->seq_header.film_grain_params_present = src->seq_header.film_grain_params_present;
    write_count += sizeof(int32_t);
    dst->seq_header.film_grain_params_present = src->seq_header.film_grain_params_present;
    write_count += sizeof(int32_t);
    dst->picture_control_set_pool_init_count = src->picture_control_set_pool_init_count;
    write_count += sizeof(int32_t);
    dst->picture_control_set_pool_init_count_child = src->picture_control_set_pool_init_count_child;
    write_count += sizeof(int32_t);
    dst->pa_reference_picture_buffer_init_count = src->pa_reference_picture_buffer_init_count;
    write_count += sizeof(int32_t);
    dst->reference_picture_buffer_init_count = src->reference_picture_buffer_init_count;
    write_count += sizeof(int32_t);
    dst->input_buffer_fifo_init_count = src->input_buffer_fifo_init_count;
    write_count += sizeof(int32_t);
    dst->overlay_input_picture_buffer_init_count = src->overlay_input_picture_buffer_init_count;
    write_count += sizeof(int32_t);
    dst->output_stream_buffer_fifo_init_count = src->output_stream_buffer_fifo_init_count;
    write_count += sizeof(int32_t);
    dst->output_recon_buffer_fifo_init_count = src->output_recon_buffer_fifo_init_count;
    write_count += sizeof(int32_t);
    dst->resource_coordination_fifo_init_count = src->resource_coordination_fifo_init_count;
    write_count += sizeof(int32_t);
    dst->picture_analysis_fifo_init_count = src->picture_analysis_fifo_init_count;
    write_count += sizeof(int32_t);
    dst->picture_decision_fifo_init_count = src->picture_decision_fifo_init_count;
    write_count += sizeof(int32_t);
    dst->motion_estimation_fifo_init_count = src->motion_estimation_fifo_init_count;
    write_count += sizeof(int32_t);
    dst->initial_rate_control_fifo_init_count = src->initial_rate_control_fifo_init_count;
    write_count += sizeof(int32_t);
    dst->picture_demux_fifo_init_count = src->picture_demux_fifo_init_count;
    write_count += sizeof(int32_t);
    dst->rate_control_tasks_fifo_init_count = src->rate_control_tasks_fifo_init_count;
    write_count += sizeof(int32_t);
    dst->rate_control_fifo_init_count = src->rate_control_fifo_init_count;
    write_count += sizeof(int32_t);
    dst->mode_decision_configuration_fifo_init_count =
        src->mode_decision_configuration_fifo_init_count;
    write_count += sizeof(int32_t);
    dst->enc_dec_fifo_init_count = src->enc_dec_fifo_init_count;
    write_count += sizeof(int32_t);
    dst->entropy_coding_fifo_init_count = src->entropy_coding_fifo_init_count;
    write_count += sizeof(int32_t);
    dst->picture_analysis_process_init_count = src->picture_analysis_process_init_count;
    write_count += sizeof(int32_t);
    dst->motion_estimation_process_init_count = src->motion_estimation_process_init_count;
    write_count += sizeof(int32_t);
    dst->source_based_operations_process_init_count =
        src->source_based_operations_process_init_count;
    write_count += sizeof(int32_t);
    dst->mode_decision_configuration_process_init_count =
        src->mode_decision_configuration_process_init_count;
    write_count += sizeof(int32_t);
    dst->enc_dec_process_init_count = src->enc_dec_process_init_count;
    write_count += sizeof(int32_t);
    dst->entropy_coding_process_init_count = src->entropy_coding_process_init_count;
    write_count += sizeof(int32_t);
    dst->total_process_init_count = src->total_process_init_count;
    write_count += sizeof(int32_t);
    dst->left_padding = src->left_padding;
    write_count += sizeof(int16_t);
    dst->right_padding = src->right_padding;
    write_count += sizeof(int16_t);
    dst->top_padding = src->top_padding;
    write_count += sizeof(int16_t);
    dst->bot_padding = src->bot_padding;
    write_count += sizeof(int16_t);
    dst->reference_count = src->reference_count;
    write_count += sizeof(uint32_t);
    for (uint8_t i = 0; i < MAX_HIERARCHICAL_LEVEL; i++) {
        dst->me_segment_column_count_array[i]   = src->me_segment_column_count_array[i];
        dst->me_segment_row_count_array[i]      = src->me_segment_row_count_array[i];
        dst->enc_dec_segment_col_count_array[i] = src->enc_dec_segment_col_count_array[i];
        dst->enc_dec_segment_row_count_array[i] = src->enc_dec_segment_row_count_array[i];
        dst->tile_group_col_count_array[i] = src->tile_group_col_count_array[i];
        dst->tile_group_row_count_array[i] = src->tile_group_row_count_array[i];
    }

    dst->cdef_segment_column_count = src->cdef_segment_column_count;
    dst->cdef_segment_row_count    = src->cdef_segment_row_count;

    dst->rest_segment_column_count      = src->rest_segment_column_count;
    dst->rest_segment_row_count         = src->rest_segment_row_count;
    dst->mrp_mode                       = src->mrp_mode;
    dst->nsq_present                    = src->nsq_present;
    dst->cdf_mode                       = src->cdf_mode;
    dst->down_sampling_method_me_search = src->down_sampling_method_me_search;
    dst->tf_segment_column_count        = src->tf_segment_column_count;
    dst->tf_segment_row_count           = src->tf_segment_row_count;
    dst->over_boundary_block_mode       = src->over_boundary_block_mode;
    dst->mfmv_enabled                   = src->mfmv_enabled;
    dst->use_input_stat_file            = src->use_input_stat_file;
    dst->use_output_stat_file           = src->use_output_stat_file;
    dst->scd_delay                      = src->scd_delay;
    return EB_ErrorNone;
}

extern EbErrorType derive_input_resolution(EbInputResolution *input_resolution, uint32_t inputSize) {
    EbErrorType return_error = EB_ErrorNone;

    if(inputSize < INPUT_SIZE_1080i_TH)
        *input_resolution = INPUT_SIZE_576p_RANGE_OR_LOWER;
    else if(inputSize < INPUT_SIZE_1080p_TH)
        *input_resolution = INPUT_SIZE_1080i_RANGE;
    else if(inputSize < INPUT_SIZE_4K_TH)
        *input_resolution = INPUT_SIZE_1080p_RANGE;
    else
        *input_resolution = INPUT_SIZE_4K_RANGE;

    return return_error;
}

static void eb_sequence_control_set_instance_dctor(EbPtr p) {
    EbSequenceControlSetInstance *obj = (EbSequenceControlSetInstance *)p;
    EB_DELETE(obj->encode_context_ptr);
    EB_DELETE(obj->scs_ptr);
    EB_DESTROY_MUTEX(obj->config_mutex);
}

EbErrorType eb_sequence_control_set_instance_ctor(EbSequenceControlSetInstance *object_ptr) {
    EbSequenceControlSetInitData scs_init_data;

    object_ptr->dctor = eb_sequence_control_set_instance_dctor;

    EB_NEW(object_ptr->encode_context_ptr, encode_context_ctor, NULL);
    scs_init_data.encode_context_ptr = object_ptr->encode_context_ptr;

    scs_init_data.sb_size = 64;

    EB_NEW(object_ptr->scs_ptr, eb_sequence_control_set_ctor, (void *)&scs_init_data);
    EB_CREATE_MUTEX(object_ptr->config_mutex);

    return EB_ErrorNone;
}

extern EbErrorType sb_params_init(SequenceControlSet *scs_ptr) {
    EbErrorType return_error = EB_ErrorNone;
    uint16_t    sb_index;
    uint16_t    raster_scan_blk_index;
    uint8_t     picture_sb_width =
        (uint8_t)((scs_ptr->seq_header.max_frame_width + scs_ptr->sb_sz - 1) / scs_ptr->sb_sz);
    uint8_t picture_sb_height =
        (uint8_t)((scs_ptr->seq_header.max_frame_height + scs_ptr->sb_sz - 1) / scs_ptr->sb_sz);
    //free old one;
    EB_FREE_ARRAY(scs_ptr->sb_params_array);

    EB_MALLOC_ARRAY(scs_ptr->sb_params_array, picture_sb_width * picture_sb_height);

    for (sb_index = 0; sb_index < picture_sb_width * picture_sb_height; ++sb_index) {
        scs_ptr->sb_params_array[sb_index].horizontal_index =
            (uint8_t)(sb_index % picture_sb_width);
        scs_ptr->sb_params_array[sb_index].vertical_index = (uint8_t)(sb_index / picture_sb_width);
        scs_ptr->sb_params_array[sb_index].origin_x =
            scs_ptr->sb_params_array[sb_index].horizontal_index * scs_ptr->sb_sz;
        scs_ptr->sb_params_array[sb_index].origin_y =
            scs_ptr->sb_params_array[sb_index].vertical_index * scs_ptr->sb_sz;

        scs_ptr->sb_params_array[sb_index].width = (uint8_t)(
            ((scs_ptr->seq_header.max_frame_width - scs_ptr->sb_params_array[sb_index].origin_x) <
             scs_ptr->sb_sz)
                ? scs_ptr->seq_header.max_frame_width - scs_ptr->sb_params_array[sb_index].origin_x
                : scs_ptr->sb_sz);

        scs_ptr->sb_params_array[sb_index].height = (uint8_t)(
            ((scs_ptr->seq_header.max_frame_height - scs_ptr->sb_params_array[sb_index].origin_y) <
             scs_ptr->sb_sz)
                ? scs_ptr->seq_header.max_frame_height - scs_ptr->sb_params_array[sb_index].origin_y
                : scs_ptr->sb_sz);

        scs_ptr->sb_params_array[sb_index].is_complete_sb =
            (uint8_t)(((scs_ptr->sb_params_array[sb_index].width == scs_ptr->sb_sz) &&
                       (scs_ptr->sb_params_array[sb_index].height == scs_ptr->sb_sz))
                          ? 1
                          : 0);

        scs_ptr->sb_params_array[sb_index].is_edge_sb =
            (scs_ptr->sb_params_array[sb_index].origin_x < scs_ptr->sb_sz) ||
                    (scs_ptr->sb_params_array[sb_index].origin_y < scs_ptr->sb_sz) ||
                    (scs_ptr->sb_params_array[sb_index].origin_x >
                     scs_ptr->seq_header.max_frame_width - scs_ptr->sb_sz) ||
                    (scs_ptr->sb_params_array[sb_index].origin_y >
                     scs_ptr->seq_header.max_frame_height - scs_ptr->sb_sz)
                ? 1
                : 0;

        for (raster_scan_blk_index = RASTER_SCAN_CU_INDEX_64x64;
             raster_scan_blk_index <= RASTER_SCAN_CU_INDEX_8x8_63;
             raster_scan_blk_index++) {
            scs_ptr->sb_params_array[sb_index].raster_scan_blk_validity[raster_scan_blk_index] =
                ((scs_ptr->sb_params_array[sb_index].origin_x +
                      raster_scan_blk_x[raster_scan_blk_index] +
                      raster_scan_blk_size[raster_scan_blk_index] >
                  scs_ptr->seq_header.max_frame_width) ||
                 (scs_ptr->sb_params_array[sb_index].origin_y +
                      raster_scan_blk_y[raster_scan_blk_index] +
                      raster_scan_blk_size[raster_scan_blk_index] >
                  scs_ptr->seq_header.max_frame_height))
                    ? EB_FALSE
                    : EB_TRUE;
        }
    }

    scs_ptr->pic_width_in_sb      = picture_sb_width;
    scs_ptr->picture_height_in_sb = picture_sb_height;
    scs_ptr->sb_total_count       = picture_sb_width * picture_sb_height;

    return return_error;
}

EbErrorType sb_geom_init(SequenceControlSet *scs_ptr) {
    uint16_t sb_index;
    uint16_t md_scan_block_index;
    uint16_t picture_sb_width =
        (scs_ptr->seq_header.max_frame_width + scs_ptr->sb_size_pix - 1) / scs_ptr->sb_size_pix;
    uint16_t picture_sb_height =
        (scs_ptr->seq_header.max_frame_height + scs_ptr->sb_size_pix - 1) / scs_ptr->sb_size_pix;

    EB_FREE_ARRAY(scs_ptr->sb_geom);
    EB_MALLOC_ARRAY(scs_ptr->sb_geom, picture_sb_width * picture_sb_height);

    for (sb_index = 0; sb_index < picture_sb_width * picture_sb_height; ++sb_index) {
        scs_ptr->sb_geom[sb_index].horizontal_index = sb_index % picture_sb_width;
        scs_ptr->sb_geom[sb_index].vertical_index   = sb_index / picture_sb_width;
        scs_ptr->sb_geom[sb_index].origin_x =
            scs_ptr->sb_geom[sb_index].horizontal_index * scs_ptr->sb_size_pix;
        scs_ptr->sb_geom[sb_index].origin_y =
            scs_ptr->sb_geom[sb_index].vertical_index * scs_ptr->sb_size_pix;

        scs_ptr->sb_geom[sb_index].width = (uint8_t)(
            ((scs_ptr->seq_header.max_frame_width - scs_ptr->sb_geom[sb_index].origin_x) <
             scs_ptr->sb_size_pix)
                ? scs_ptr->seq_header.max_frame_width - scs_ptr->sb_geom[sb_index].origin_x
                : scs_ptr->sb_size_pix);

        scs_ptr->sb_geom[sb_index].height = (uint8_t)(
            ((scs_ptr->seq_header.max_frame_height - scs_ptr->sb_geom[sb_index].origin_y) <
             scs_ptr->sb_size_pix)
                ? scs_ptr->seq_header.max_frame_height - scs_ptr->sb_geom[sb_index].origin_y
                : scs_ptr->sb_size_pix);

        scs_ptr->sb_geom[sb_index].is_complete_sb =
            (uint8_t)(((scs_ptr->sb_geom[sb_index].width == scs_ptr->sb_size_pix) &&
                       (scs_ptr->sb_geom[sb_index].height == scs_ptr->sb_size_pix))
                          ? 1
                          : 0);

        uint16_t max_block_count = scs_ptr->max_block_cnt;

        for (md_scan_block_index = 0; md_scan_block_index < max_block_count;
             md_scan_block_index++) {
            const BlockGeom *blk_geom = get_blk_geom_mds(md_scan_block_index);
            if (scs_ptr->over_boundary_block_mode == 1) {
                scs_ptr->sb_geom[sb_index].block_is_allowed[md_scan_block_index] =
                    ((scs_ptr->sb_geom[sb_index].origin_x + blk_geom->origin_x +
                          blk_geom->bwidth / 2 <
                      scs_ptr->seq_header.max_frame_width) &&
                     (scs_ptr->sb_geom[sb_index].origin_y + blk_geom->origin_y +
                          blk_geom->bheight / 2 <
                      scs_ptr->seq_header.max_frame_height))
                        ? EB_TRUE
                        : EB_FALSE;

                // Temporary if the cropped width is not 4, 8, 16, 32, 64 and 128, the block is not allowed. To be removed after intrinsic functions for NxM spatial_full_distortion_kernel_func_ptr_array are added
                int32_t cropped_width =
                    MIN(blk_geom->bwidth,
                        scs_ptr->seq_header.max_frame_width -
                            (scs_ptr->sb_geom[sb_index].origin_x + blk_geom->origin_x));
                if (cropped_width != 4 && cropped_width != 8 && cropped_width != 16 &&
                    cropped_width != 32 && cropped_width != 64 && cropped_width != 128)
                    scs_ptr->sb_geom[sb_index].block_is_allowed[md_scan_block_index] = EB_FALSE;

                if (blk_geom->shape != PART_N) blk_geom = get_blk_geom_mds(blk_geom->sqi_mds);
                scs_ptr->sb_geom[sb_index].block_is_inside_md_scan[md_scan_block_index] =
                    ((scs_ptr->sb_geom[sb_index].origin_x >= scs_ptr->seq_header.max_frame_width) ||
                     (scs_ptr->sb_geom[sb_index].origin_y >= scs_ptr->seq_header.max_frame_height))
                        ? EB_FALSE
                        : EB_TRUE;
            } else {
                if (blk_geom->shape != PART_N) blk_geom = get_blk_geom_mds(blk_geom->sqi_mds);

                scs_ptr->sb_geom[sb_index].block_is_allowed[md_scan_block_index] =
                    ((scs_ptr->sb_geom[sb_index].origin_x + blk_geom->origin_x + blk_geom->bwidth >
                      scs_ptr->seq_header.max_frame_width) ||
                     (scs_ptr->sb_geom[sb_index].origin_y + blk_geom->origin_y + blk_geom->bheight >
                      scs_ptr->seq_header.max_frame_height))
                        ? EB_FALSE
                        : EB_TRUE;

                scs_ptr->sb_geom[sb_index].block_is_inside_md_scan[md_scan_block_index] =
                    ((scs_ptr->sb_geom[sb_index].origin_x + blk_geom->origin_x + blk_geom->bwidth >
                      scs_ptr->seq_header.max_frame_width) ||
                     (scs_ptr->sb_geom[sb_index].origin_y + blk_geom->origin_y + blk_geom->bheight >
                      scs_ptr->seq_header.max_frame_height))
                        ? EB_FALSE
                        : EB_TRUE;
            }
        }
    }

    scs_ptr->sb_tot_cnt = picture_sb_width * picture_sb_height;

    return 0;
}
