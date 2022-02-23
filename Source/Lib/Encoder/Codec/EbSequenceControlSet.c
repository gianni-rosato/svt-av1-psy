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

#include "EbSequenceControlSet.h"
#include "EbUtility.h"

static void svt_sequence_control_set_dctor(EbPtr p) {
    SequenceControlSet *obj = (SequenceControlSet *)p;
    EB_FREE_ARRAY(obj->sb_params_array);
    EB_FREE_ARRAY(obj->sb_geom);
}

/**************************************************************************************************
    General notes on how Sequence Control Sets (SCS) are used.

    SequenceControlSetInstance
        is the primary copy that interacts with the API in real-time.  When a
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
EbErrorType svt_sequence_control_set_ctor(SequenceControlSet *scs_ptr, EbPtr object_init_data_ptr) {
    EbSequenceControlSetInitData *scs_init_data = (EbSequenceControlSetInitData *)
        object_init_data_ptr;
    uint32_t segment_index;
    scs_ptr->mvrate_set = 0;
    scs_ptr->dctor      = svt_sequence_control_set_dctor;

    scs_ptr->sb_sz = 64;

    // Segments
    for (segment_index = 0; segment_index < MAX_TEMPORAL_LAYERS; ++segment_index) {
        scs_ptr->me_segment_column_count_array[segment_index]   = 1;
        scs_ptr->me_segment_row_count_array[segment_index]      = 1;
        scs_ptr->enc_dec_segment_col_count_array[segment_index] = 1;
        scs_ptr->enc_dec_segment_row_count_array[segment_index] = 1;
    }
    scs_ptr->tpl_segment_col_count_array = 1;
    scs_ptr->tpl_segment_row_count_array = 1;

    // Encode Context
    if (scs_init_data != NULL)
        scs_ptr->encode_context_ptr = scs_init_data->encode_context_ptr;

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
    scs_ptr->sb_sz            = 64;
    scs_ptr->super_block_size = 128;
    scs_ptr->max_sb_depth     = 3;

    // CU
    scs_ptr->max_blk_size   = 64;
    scs_ptr->min_blk_size   = 8;
    scs_ptr->max_intra_size = 32;
    scs_ptr->min_intra_size = 8;
    // Quantization
    scs_ptr->static_config.qp = 20;
    // AV1 features
    scs_ptr->enable_warped_motion        = DEFAULT;
    scs_ptr->enable_global_motion        = 1;
    scs_ptr->sg_filter_mode              = DEFAULT;
    scs_ptr->wn_filter_mode              = DEFAULT;
    scs_ptr->inter_intra_compound        = DEFAULT;
    scs_ptr->enable_paeth                = DEFAULT;
    scs_ptr->enable_smooth               = DEFAULT;
    scs_ptr->spatial_sse_full_loop_level = DEFAULT;
    scs_ptr->over_bndry_blk              = DEFAULT;
    scs_ptr->new_nearest_comb_inject     = DEFAULT;
    scs_ptr->nsq_table                   = DEFAULT;
    scs_ptr->frame_end_cdf_update        = DEFAULT;
    scs_ptr->set_chroma_mode             = DEFAULT;
    scs_ptr->disable_cfl_flag            = DEFAULT;
    scs_ptr->obmc_level                  = DEFAULT;
    scs_ptr->rdoq_level                  = DEFAULT;
    scs_ptr->pred_me                     = DEFAULT;
    scs_ptr->bipred_3x3_inject           = DEFAULT;
    scs_ptr->compound_level              = DEFAULT;
    scs_ptr->filter_intra_level          = DEFAULT;
    scs_ptr->enable_intra_edge_filter    = DEFAULT;
    scs_ptr->pic_based_rate_est          = DEFAULT;
    scs_ptr->ext_block_flag              = EB_FALSE;
    scs_ptr->ten_bit_format              = EB_FALSE;

    // Initialize SB params
    //allocation will happen in ress-corrd
    scs_ptr->sb_params_array = 0;

    scs_ptr->seq_header.frame_width_bits              = 16;
    scs_ptr->seq_header.frame_height_bits             = 16;
    scs_ptr->seq_header.frame_id_numbers_present_flag = 0;
    //    cm->large_scale_tile ? 0 : cm->error_resilient_mode;
    scs_ptr->seq_header.frame_id_length       = FRAME_ID_LENGTH;
    scs_ptr->seq_header.delta_frame_id_length = DELTA_FRAME_ID_LENGTH;

    if (scs_init_data && scs_init_data->sb_size == 128) {
        scs_ptr->seq_header.sb_size      = BLOCK_128X128;
        scs_ptr->sb_size_pix             = 128;
        scs_ptr->seq_header.sb_mi_size   = 32; // Size of the superblock in units of MI blocks
        scs_ptr->seq_header.sb_size_log2 = 5;
    } else {
        scs_ptr->seq_header.sb_size = BLOCK_64X64;
        scs_ptr->sb_size_pix        = 64;

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

    scs_ptr->seq_header.order_hint_info.order_hint_bits = 7;

    scs_ptr->seq_header.seq_force_screen_content_tools = 2;
    // 0 - force off
    // 1 - force on
    // 2 - adaptive
    scs_ptr->seq_header.seq_force_integer_mv = 2; // 0 - Not to force. MV can be in 1/4 or 1/8
    // 1 - force to integer
    // 2 - adaptive

    scs_ptr->seq_header.order_hint_info.enable_ref_frame_mvs = 1;
    if (scs_ptr->static_config.cdef_level == DEFAULT)
        scs_ptr->seq_header.cdef_level = 1;
    else
        scs_ptr->seq_header.cdef_level = (uint8_t)scs_ptr->static_config.cdef_level;

    if (scs_ptr->static_config.enable_restoration_filtering == DEFAULT)
        scs_ptr->seq_header.enable_restoration = 1;
    else
        scs_ptr->seq_header.enable_restoration =
            (uint8_t)scs_ptr->static_config.enable_restoration_filtering;

    if (scs_ptr->enable_intra_edge_filter == DEFAULT)
        scs_ptr->seq_header.enable_intra_edge_filter = 1;
    else
        scs_ptr->seq_header.enable_intra_edge_filter = (uint8_t)scs_ptr->enable_intra_edge_filter;

    if (scs_ptr->pic_based_rate_est == DEFAULT)
        scs_ptr->seq_header.pic_based_rate_est = 0;
    else
        scs_ptr->seq_header.pic_based_rate_est = (uint8_t)scs_ptr->pic_based_rate_est;
    if (scs_ptr->enable_warped_motion == DEFAULT)
        scs_ptr->seq_header.enable_warped_motion = 1;
    else
        scs_ptr->seq_header.enable_warped_motion = (uint8_t)scs_ptr->enable_warped_motion;

    scs_ptr->film_grain_random_seed = 7391;
    scs_ptr->reference_count        = 4;

    // Set the block mean calculation prec
    scs_ptr->block_mean_calc_prec = BLOCK_MEAN_PREC_SUB;

    // Set Picture Parameters for statistics gathering
    scs_ptr->picture_analysis_number_of_regions_per_width =
        HIGHER_THAN_CLASS_1_REGION_SPLIT_PER_WIDTH;
    scs_ptr->picture_analysis_number_of_regions_per_height =
        HIGHER_THAN_CLASS_1_REGION_SPLIT_PER_HEIGHT;

    scs_ptr->palette_level = DEFAULT;

    // intra angle delta
    scs_ptr->intra_angle_delta = DEFAULT;

    scs_ptr->is_16bit_pipeline = EB_FALSE;

    scs_ptr->intrabc_mode = DEFAULT;

    scs_ptr->enable_qp_scaling_flag = 1;

    scs_ptr->enable_adaptive_mini_gop = 0;

    scs_ptr->max_heirachical_level = 5;

    scs_ptr->speed_control_flag = 0;

    return EB_ErrorNone;
}

EbErrorType svt_sequence_control_set_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr) {
    SequenceControlSet *obj;

    *object_dbl_ptr = NULL;
    EB_NEW(obj, svt_sequence_control_set_ctor, object_init_data_ptr);
    *object_dbl_ptr = obj;
    return EB_ErrorNone;
}

/************************************************
 * Sequence Control Set Copy
 ************************************************/
EbErrorType copy_sequence_control_set(SequenceControlSet *dst, SequenceControlSet *src) {
    dst->static_config                = src->static_config;
    dst->encode_context_ptr           = src->encode_context_ptr;
    dst->chroma_format_idc            = src->chroma_format_idc;
    dst->max_temporal_layers          = src->max_temporal_layers;
    dst->bits_for_picture_order_count = src->bits_for_picture_order_count;
    dst->max_input_luma_width         = src->max_input_luma_width;
    dst->max_input_luma_height        = src->max_input_luma_height;
    dst->max_input_chroma_height      = src->max_input_chroma_height;
    dst->max_input_chroma_width       = src->max_input_chroma_width;
    dst->max_input_pad_right          = src->max_input_pad_right;
    dst->max_input_pad_bottom         = src->max_input_pad_bottom;
    dst->seq_header.max_frame_width   = src->seq_header.max_frame_width;
    dst->seq_header.max_frame_height  = src->seq_header.max_frame_height;
    dst->chroma_width                 = src->chroma_width;
    dst->chroma_height                = src->chroma_height;
    dst->pad_right                    = src->pad_right;
    dst->pad_bottom                   = src->pad_bottom;
    dst->frame_rate                   = src->frame_rate;
    //dst->input_bitdepth = src->input_bitdepth;
    //dst->output_bitdepth = src->output_bitdepth;
    dst->encoder_bit_depth                         = src->encoder_bit_depth;
    dst->subsampling_x                             = src->subsampling_x;
    dst->subsampling_y                             = src->subsampling_y;
    dst->pred_struct_ptr                           = src->pred_struct_ptr;
    dst->max_ref_count                             = src->max_ref_count;
    dst->super_block_size                          = src->super_block_size;
    dst->sb_sz                                     = src->sb_sz;
    dst->max_sb_depth                              = src->max_sb_depth;
    dst->max_blk_size                              = src->max_blk_size;
    dst->min_blk_size                              = src->min_blk_size;
    dst->max_intra_size                            = src->max_intra_size;
    dst->min_intra_size                            = src->min_intra_size;
    dst->static_config.qp                          = src->static_config.qp;
    dst->seq_header.film_grain_params_present      = src->seq_header.film_grain_params_present;
    dst->picture_control_set_pool_init_count       = src->picture_control_set_pool_init_count;
    dst->me_pool_init_count                        = src->me_pool_init_count;
    dst->picture_control_set_pool_init_count_child = src->picture_control_set_pool_init_count_child;
    dst->enc_dec_pool_init_count                   = src->enc_dec_pool_init_count;
    dst->pa_reference_picture_buffer_init_count    = src->pa_reference_picture_buffer_init_count;
    dst->reference_picture_buffer_init_count       = src->reference_picture_buffer_init_count;
    dst->input_buffer_fifo_init_count              = src->input_buffer_fifo_init_count;
    dst->overlay_input_picture_buffer_init_count   = src->overlay_input_picture_buffer_init_count;
    dst->output_stream_buffer_fifo_init_count      = src->output_stream_buffer_fifo_init_count;
    dst->output_recon_buffer_fifo_init_count       = src->output_recon_buffer_fifo_init_count;
    dst->resource_coordination_fifo_init_count     = src->resource_coordination_fifo_init_count;
    dst->picture_analysis_fifo_init_count          = src->picture_analysis_fifo_init_count;
    dst->picture_decision_fifo_init_count          = src->picture_decision_fifo_init_count;
    dst->motion_estimation_fifo_init_count         = src->motion_estimation_fifo_init_count;
    dst->initial_rate_control_fifo_init_count      = src->initial_rate_control_fifo_init_count;
    dst->picture_demux_fifo_init_count             = src->picture_demux_fifo_init_count;
    dst->rate_control_tasks_fifo_init_count        = src->rate_control_tasks_fifo_init_count;
    dst->rate_control_fifo_init_count              = src->rate_control_fifo_init_count;
    dst->mode_decision_configuration_fifo_init_count =
        src->mode_decision_configuration_fifo_init_count;
    dst->enc_dec_fifo_init_count              = src->enc_dec_fifo_init_count;
    dst->entropy_coding_fifo_init_count       = src->entropy_coding_fifo_init_count;
    dst->picture_analysis_process_init_count  = src->picture_analysis_process_init_count;
    dst->motion_estimation_process_init_count = src->motion_estimation_process_init_count;
    dst->source_based_operations_process_init_count =
        src->source_based_operations_process_init_count;
    dst->mode_decision_configuration_process_init_count =
        src->mode_decision_configuration_process_init_count;
    dst->enc_dec_process_init_count        = src->enc_dec_process_init_count;
    dst->entropy_coding_process_init_count = src->entropy_coding_process_init_count;
    dst->total_process_init_count          = src->total_process_init_count;
    dst->left_padding                      = src->left_padding;
    dst->right_padding                     = src->right_padding;
    dst->top_padding                       = src->top_padding;
    dst->bot_padding                       = src->bot_padding;
    dst->reference_count                   = src->reference_count;
    for (uint8_t i = 0; i < MAX_HIERARCHICAL_LEVEL; i++) {
        dst->me_segment_column_count_array[i]   = src->me_segment_column_count_array[i];
        dst->me_segment_row_count_array[i]      = src->me_segment_row_count_array[i];
        dst->enc_dec_segment_col_count_array[i] = src->enc_dec_segment_col_count_array[i];
        dst->enc_dec_segment_row_count_array[i] = src->enc_dec_segment_row_count_array[i];
        dst->tile_group_col_count_array[i]      = src->tile_group_col_count_array[i];
        dst->tile_group_row_count_array[i]      = src->tile_group_row_count_array[i];
    }

    dst->tpl_segment_col_count_array = src->tpl_segment_col_count_array;
    dst->tpl_segment_row_count_array = src->tpl_segment_row_count_array;

    dst->cdef_segment_column_count = src->cdef_segment_column_count;
    dst->cdef_segment_row_count    = src->cdef_segment_row_count;

    dst->rest_segment_column_count      = src->rest_segment_column_count;
    dst->rest_segment_row_count         = src->rest_segment_row_count;
    dst->cdf_mode                       = src->cdf_mode;
    dst->down_sampling_method_me_search = src->down_sampling_method_me_search;
    dst->tf_segment_column_count        = src->tf_segment_column_count;
    dst->tf_segment_row_count           = src->tf_segment_row_count;
    dst->over_boundary_block_mode       = src->over_boundary_block_mode;
    dst->mfmv_enabled                   = src->mfmv_enabled;
    dst->scd_delay                      = src->scd_delay;
    dst->in_loop_ois                    = src->in_loop_ois;
    dst->enable_pic_mgr_dec_order       = src->enable_pic_mgr_dec_order;
    dst->enable_dec_order               = src->enable_dec_order;
    dst->lap_enabled                    = src->lap_enabled;

    dst->lad_mg                        = src->lad_mg;
    dst->tpl_lad_mg                    = src->tpl_lad_mg;
    dst->use_boundaries_in_rest_search = src->use_boundaries_in_rest_search;
    dst->geom_idx                      = src->geom_idx;
    dst->max_block_cnt                 = src->max_block_cnt;
    dst->rc_stat_gen_pass_mode         = src->rc_stat_gen_pass_mode;
    dst->cqp_base_q_tf                 = src->cqp_base_q_tf;
    dst->cqp_base_q                    = src->cqp_base_q;
    dst->is_short_clip                 = src->is_short_clip;
    dst->input_resolution              = src->input_resolution;
    dst->tf_params_per_type[0]         = src->tf_params_per_type[0];
    dst->tf_params_per_type[1]         = src->tf_params_per_type[1];
    dst->tf_params_per_type[2]         = src->tf_params_per_type[2];
    dst->vq_ctrls                      = src->vq_ctrls;
    dst->mrp_ctrls                     = src->mrp_ctrls;
    dst->passes                        = src->passes;
    dst->ipp_pass_ctrls                = src->ipp_pass_ctrls;
    dst->mid_pass_ctrls                = src->mid_pass_ctrls;
    dst->ipp_was_ds                    = src->ipp_was_ds;
    dst->final_pass_preset             = src->final_pass_preset;
    dst->palette_level                 = src->palette_level;
    dst->intra_angle_delta             = src->intra_angle_delta;
    dst->is_16bit_pipeline             = src->is_16bit_pipeline;
    dst->enable_warped_motion          = src->enable_warped_motion;
    dst->enable_global_motion          = src->enable_global_motion;
    dst->sg_filter_mode                = src->sg_filter_mode;
    dst->wn_filter_mode                = src->wn_filter_mode;
    dst->inter_intra_compound          = src->inter_intra_compound;
    dst->enable_paeth                  = src->enable_paeth;
    dst->enable_smooth                 = src->enable_smooth;
    dst->spatial_sse_full_loop_level   = src->spatial_sse_full_loop_level;
    dst->over_bndry_blk                = src->spatial_sse_full_loop_level;
    dst->new_nearest_comb_inject       = src->new_nearest_comb_inject;
    dst->nsq_table                     = src->nsq_table;
    dst->frame_end_cdf_update          = src->frame_end_cdf_update;
    dst->set_chroma_mode               = src->set_chroma_mode;
    dst->disable_cfl_flag              = src->disable_cfl_flag;
    dst->obmc_level                    = src->obmc_level;
    dst->rdoq_level                    = src->rdoq_level;
    dst->pred_me                       = src->pred_me;
    dst->bipred_3x3_inject             = src->bipred_3x3_inject;
    dst->compound_level                = src->compound_level;
    dst->filter_intra_level            = src->filter_intra_level;
    dst->enable_intra_edge_filter      = src->enable_intra_edge_filter;
    dst->pic_based_rate_est            = src->pic_based_rate_est;
    dst->ext_block_flag                = src->ext_block_flag;
    dst->intrabc_mode                  = src->intrabc_mode;
    dst->enable_hbd_mode_decision      = src->enable_hbd_mode_decision;
    dst->enable_qp_scaling_flag        = src->enable_qp_scaling_flag;
    dst->ten_bit_format                = src->ten_bit_format;
    dst->enable_adaptive_mini_gop      = src->enable_adaptive_mini_gop;
    dst->max_heirachical_level         = src->max_heirachical_level;
    dst->speed_control_flag            = src->speed_control_flag;
    dst->tpl_level                     = src->tpl_level;
    dst->calculate_variance            = src->calculate_variance;
    return EB_ErrorNone;
}

extern EbErrorType derive_input_resolution(EbInputResolution *input_resolution,
                                           uint32_t           inputSize) {
    EbErrorType return_error = EB_ErrorNone;
    if (inputSize < INPUT_SIZE_240p_TH)
        *input_resolution = INPUT_SIZE_240p_RANGE;
    else if (inputSize < INPUT_SIZE_360p_TH)
        *input_resolution = INPUT_SIZE_360p_RANGE;
    else if (inputSize < INPUT_SIZE_480p_TH)
        *input_resolution = INPUT_SIZE_480p_RANGE;
    else if (inputSize < INPUT_SIZE_720p_TH)
        *input_resolution = INPUT_SIZE_720p_RANGE;
    else if (inputSize < INPUT_SIZE_1080p_TH)
        *input_resolution = INPUT_SIZE_1080p_RANGE;
    else if (inputSize < INPUT_SIZE_4K_TH)
        *input_resolution = INPUT_SIZE_4K_RANGE;
    else
        *input_resolution = INPUT_SIZE_8K_RANGE;

    return return_error;
}

static void svt_sequence_control_set_instance_dctor(EbPtr p) {
    EbSequenceControlSetInstance *obj = (EbSequenceControlSetInstance *)p;
    EB_DELETE(obj->encode_context_ptr);
    EB_DELETE(obj->scs_ptr);
    EB_DESTROY_MUTEX(obj->config_mutex);
}

EbErrorType svt_sequence_control_set_instance_ctor(EbSequenceControlSetInstance *object_ptr) {
    EbSequenceControlSetInitData scs_init_data;

    object_ptr->dctor = svt_sequence_control_set_instance_dctor;

    EB_NEW(object_ptr->encode_context_ptr, encode_context_ctor, NULL);
    scs_init_data.encode_context_ptr = object_ptr->encode_context_ptr;

    scs_init_data.sb_size = 64;

    EB_NEW(object_ptr->scs_ptr, svt_sequence_control_set_ctor, (void *)&scs_init_data);
    EB_CREATE_MUTEX(object_ptr->config_mutex);

    return EB_ErrorNone;
}

extern EbErrorType sb_params_init(SequenceControlSet *scs_ptr) {
    EbErrorType return_error = EB_ErrorNone;
    uint16_t    sb_index;
    uint16_t    raster_scan_blk_index;

    uint16_t picture_sb_width = (scs_ptr->seq_header.max_frame_width + scs_ptr->sb_sz - 1) /
        scs_ptr->sb_sz;
    uint16_t picture_sb_height = (scs_ptr->seq_header.max_frame_height + scs_ptr->sb_sz - 1) /
        scs_ptr->sb_sz;
    //free old one;
    EB_FREE_ARRAY(scs_ptr->sb_params_array);

    EB_MALLOC_ARRAY(scs_ptr->sb_params_array, picture_sb_width * picture_sb_height);

    for (sb_index = 0; sb_index < picture_sb_width * picture_sb_height; ++sb_index) {
        scs_ptr->sb_params_array[sb_index].horizontal_index = (uint8_t)(sb_index %
                                                                        picture_sb_width);
        scs_ptr->sb_params_array[sb_index].vertical_index = (uint8_t)(sb_index / picture_sb_width);
        scs_ptr->sb_params_array[sb_index].origin_x =
            scs_ptr->sb_params_array[sb_index].horizontal_index * scs_ptr->sb_sz;
        scs_ptr->sb_params_array[sb_index].origin_y =
            scs_ptr->sb_params_array[sb_index].vertical_index * scs_ptr->sb_sz;

        scs_ptr->sb_params_array[sb_index].width =
            (uint8_t)(((scs_ptr->seq_header.max_frame_width -
                        scs_ptr->sb_params_array[sb_index].origin_x) < scs_ptr->sb_sz)
                          ? scs_ptr->seq_header.max_frame_width -
                              scs_ptr->sb_params_array[sb_index].origin_x
                          : scs_ptr->sb_sz);

        scs_ptr->sb_params_array[sb_index].height =
            (uint8_t)(((scs_ptr->seq_header.max_frame_height -
                        scs_ptr->sb_params_array[sb_index].origin_y) < scs_ptr->sb_sz)
                          ? scs_ptr->seq_header.max_frame_height -
                              scs_ptr->sb_params_array[sb_index].origin_y
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

EbErrorType rtime_alloc_sb_geom(SequenceControlSet *scs_ptr, uint32_t size) {
    EB_MALLOC_ARRAY(scs_ptr->sb_geom, size);
    return EB_ErrorNone;
}
EbErrorType sb_geom_init(SequenceControlSet *scs_ptr) {
    uint16_t sb_index;
    uint16_t md_scan_block_index;
    uint16_t picture_sb_width = (scs_ptr->seq_header.max_frame_width + scs_ptr->sb_size_pix - 1) /
        scs_ptr->sb_size_pix;
    uint16_t picture_sb_height = (scs_ptr->seq_header.max_frame_height + scs_ptr->sb_size_pix - 1) /
        scs_ptr->sb_size_pix;

    EB_FREE_ARRAY(scs_ptr->sb_geom);
    rtime_alloc_sb_geom(scs_ptr, picture_sb_width * picture_sb_height);

    for (sb_index = 0; sb_index < picture_sb_width * picture_sb_height; ++sb_index) {
        scs_ptr->sb_geom[sb_index].horizontal_index = sb_index % picture_sb_width;
        scs_ptr->sb_geom[sb_index].vertical_index   = sb_index / picture_sb_width;
        scs_ptr->sb_geom[sb_index].origin_x         = scs_ptr->sb_geom[sb_index].horizontal_index *
            scs_ptr->sb_size_pix;
        scs_ptr->sb_geom[sb_index].origin_y = scs_ptr->sb_geom[sb_index].vertical_index *
            scs_ptr->sb_size_pix;

        scs_ptr->sb_geom[sb_index].width = (uint8_t)(((scs_ptr->seq_header.max_frame_width -
                                                       scs_ptr->sb_geom[sb_index].origin_x) <
                                                      scs_ptr->sb_size_pix)
                                                         ? scs_ptr->seq_header.max_frame_width -
                                                             scs_ptr->sb_geom[sb_index].origin_x
                                                         : scs_ptr->sb_size_pix);

        scs_ptr->sb_geom[sb_index].height = (uint8_t)(((scs_ptr->seq_header.max_frame_height -
                                                        scs_ptr->sb_geom[sb_index].origin_y) <
                                                       scs_ptr->sb_size_pix)
                                                          ? scs_ptr->seq_header.max_frame_height -
                                                              scs_ptr->sb_geom[sb_index].origin_y
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

                if (blk_geom->shape != PART_N)
                    blk_geom = get_blk_geom_mds(blk_geom->sqi_mds);
                scs_ptr->sb_geom[sb_index].block_is_inside_md_scan[md_scan_block_index] =
                    ((scs_ptr->sb_geom[sb_index].origin_x >= scs_ptr->seq_header.max_frame_width) ||
                     (scs_ptr->sb_geom[sb_index].origin_y >= scs_ptr->seq_header.max_frame_height))
                    ? EB_FALSE
                    : EB_TRUE;
            } else {
                if (blk_geom->shape != PART_N)
                    blk_geom = get_blk_geom_mds(blk_geom->sqi_mds);

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
