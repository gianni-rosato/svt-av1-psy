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

/***************************************
 * Includes
 ***************************************/

#include <stdlib.h>
#include <string.h>

#include "EbAppContext.h"
#include "EbAppConfig.h"

#define IS_16_BIT(bit_depth) (bit_depth == 10 ? 1 : 0)

/***************************************
 * Variables Defining a memory table
 *  hosting all allocated pointers
 ***************************************/
EbMemoryMapEntry *       app_memory_map;
uint32_t *               app_memory_map_index;
uint64_t *               total_app_memory;
uint32_t                 app_malloc_count = 0;
static EbMemoryMapEntry *app_memory_map_all_channels[MAX_CHANNEL_NUMBER];
static uint32_t          app_memory_map_index_all_channels[MAX_CHANNEL_NUMBER];
static uint64_t          app_memory_mallocd_all_channels[MAX_CHANNEL_NUMBER];

/***************************************
* Allocation and initializing a memory table
*  hosting all allocated pointers
***************************************/
void allocate_memory_table(uint32_t instance_idx) {
    // Malloc Memory Table for the instance @ instance_idx
    app_memory_map_all_channels[instance_idx] =
        (EbMemoryMapEntry *)malloc(sizeof(EbMemoryMapEntry) * MAX_APP_NUM_PTR);

    // Init the table index
    app_memory_map_index_all_channels[instance_idx] = 0;

    // Size of the table
    app_memory_mallocd_all_channels[instance_idx] = sizeof(EbMemoryMapEntry) * MAX_APP_NUM_PTR;
    total_app_memory                              = &app_memory_mallocd_all_channels[instance_idx];

    // Set pointer to the first entry
    app_memory_map = app_memory_map_all_channels[instance_idx];

    // Set index to the first entry
    app_memory_map_index = &app_memory_map_index_all_channels[instance_idx];

    // Init Number of pointers
    app_malloc_count = 0;

    return;
}

/*************************************
**************************************
*** Helper functions Input / Output **
**************************************
**************************************/

/***********************************************
* Copy configuration parameters from
*  The config structure, to the
*  callback structure to send to the library
***********************************************/
EbErrorType copy_configuration_parameters(EbConfig *config, EbAppContext *callback_data,
                                          uint32_t instance_idx) {
    EbErrorType return_error = EB_ErrorNone;
    uint32_t    hme_region_index;

    // Assign Instance index to the library
    callback_data->instance_idx = (uint8_t)instance_idx;
    // Initialize Port Activity Flags
    callback_data->output_stream_port_active                = APP_PortActive;
    callback_data->eb_enc_parameters.source_width           = config->source_width;
    callback_data->eb_enc_parameters.source_height          = config->source_height;
    callback_data->eb_enc_parameters.render_width           = config->input_padded_width;
    callback_data->eb_enc_parameters.render_height          = config->input_padded_height;
    callback_data->eb_enc_parameters.intra_period_length    = config->intra_period;
    callback_data->eb_enc_parameters.intra_refresh_type     = config->intra_refresh_type;
    callback_data->eb_enc_parameters.enc_mode               = (EbBool)config->enc_mode;
    callback_data->eb_enc_parameters.frame_rate             = config->frame_rate;
    callback_data->eb_enc_parameters.frame_rate_denominator = config->frame_rate_denominator;
    callback_data->eb_enc_parameters.frame_rate_numerator   = config->frame_rate_numerator;
    callback_data->eb_enc_parameters.hierarchical_levels    = config->hierarchical_levels;
    callback_data->eb_enc_parameters.pred_structure         = (uint8_t)config->pred_structure;
    callback_data->eb_enc_parameters.ext_block_flag         = config->ext_block_flag;
    callback_data->eb_enc_parameters.tile_rows              = config->tile_rows;
    callback_data->eb_enc_parameters.tile_columns           = config->tile_columns;
    callback_data->eb_enc_parameters.scene_change_detection = config->scene_change_detection;
    callback_data->eb_enc_parameters.look_ahead_distance    = config->look_ahead_distance;
    callback_data->eb_enc_parameters.enable_tpl_la          = config->enable_tpl_la;
    callback_data->eb_enc_parameters.rate_control_mode      = config->rate_control_mode;
    callback_data->eb_enc_parameters.target_bit_rate        = config->target_bit_rate;
    callback_data->eb_enc_parameters.max_qp_allowed         = config->max_qp_allowed;
    callback_data->eb_enc_parameters.min_qp_allowed         = config->min_qp_allowed;
    callback_data->eb_enc_parameters.vbv_bufsize            = config->vbv_bufsize;
    callback_data->eb_enc_parameters.vbr_bias_pct           = config->vbr_bias_pct;
    callback_data->eb_enc_parameters.vbr_min_section_pct    = config->vbr_min_section_pct;
    callback_data->eb_enc_parameters.vbr_max_section_pct    = config->vbr_max_section_pct;
    callback_data->eb_enc_parameters.under_shoot_pct        = config->under_shoot_pct;
    callback_data->eb_enc_parameters.over_shoot_pct         = config->over_shoot_pct;
    callback_data->eb_enc_parameters.enable_adaptive_quantization =
        (EbBool)config->enable_adaptive_quantization;
    callback_data->eb_enc_parameters.qp                   = config->qp;
    callback_data->eb_enc_parameters.use_qp_file          = (EbBool)config->use_qp_file;
    callback_data->eb_enc_parameters.rc_twopass_stats_in  = config->rc_twopass_stats_in;
    callback_data->eb_enc_parameters.rc_firstpass_stats_out = config->rc_firstpass_stats_out;
    callback_data->eb_enc_parameters.stat_report          = (EbBool)config->stat_report;
    callback_data->eb_enc_parameters.disable_dlf_flag     = (EbBool)config->disable_dlf_flag;
    callback_data->eb_enc_parameters.enable_warped_motion = config->enable_warped_motion;
    callback_data->eb_enc_parameters.enable_global_motion = (EbBool)config->enable_global_motion;
    callback_data->eb_enc_parameters.cdef_level               = config->cdef_level;
    callback_data->eb_enc_parameters.enable_restoration_filtering = config->enable_restoration_filtering;
    callback_data->eb_enc_parameters.sg_filter_mode           = config->sg_filter_mode;
    callback_data->eb_enc_parameters.wn_filter_mode           = config->wn_filter_mode;
    callback_data->eb_enc_parameters.intra_angle_delta        = config->intra_angle_delta;
    callback_data->eb_enc_parameters.inter_intra_compound     = config->inter_intra_compound;
    callback_data->eb_enc_parameters.enable_paeth             = config->enable_paeth;
    callback_data->eb_enc_parameters.enable_smooth            = config->enable_smooth;
    callback_data->eb_enc_parameters.mrp_level                = config->mrp_level;
    callback_data->eb_enc_parameters.enable_mfmv              = config->enable_mfmv;
    callback_data->eb_enc_parameters.enable_redundant_blk     = config->enable_redundant_blk;
    callback_data->eb_enc_parameters.spatial_sse_full_loop_level = config->spatial_sse_full_loop_level;
    callback_data->eb_enc_parameters.over_bndry_blk           = config->over_bndry_blk;
    callback_data->eb_enc_parameters.new_nearest_comb_inject  = config->new_nearest_comb_inject;
    callback_data->eb_enc_parameters.nsq_table                = config->nsq_table;
    callback_data->eb_enc_parameters.frame_end_cdf_update     = config->frame_end_cdf_update;
    callback_data->eb_enc_parameters.obmc_level               = (int8_t)config->obmc_level;
    callback_data->eb_enc_parameters.rdoq_level               = config->rdoq_level;
    callback_data->eb_enc_parameters.pred_me                  = config->pred_me;
    callback_data->eb_enc_parameters.bipred_3x3_inject        = config->bipred_3x3_inject;
    callback_data->eb_enc_parameters.compound_level           = config->compound_level;
    callback_data->eb_enc_parameters.set_chroma_mode          = config->set_chroma_mode;
    callback_data->eb_enc_parameters.disable_cfl_flag         = config->disable_cfl_flag;
    callback_data->eb_enc_parameters.filter_intra_level       = (int8_t)config->filter_intra_level;
    callback_data->eb_enc_parameters.pic_based_rate_est       = config->pic_based_rate_est;
    callback_data->eb_enc_parameters.enable_intra_edge_filter = config->enable_intra_edge_filter;
    callback_data->eb_enc_parameters.use_default_me_hme       = (EbBool)config->use_default_me_hme;
    callback_data->eb_enc_parameters.enable_hme_flag          = (EbBool)config->enable_hme_flag;
    callback_data->eb_enc_parameters.enable_hme_level0_flag =
        (EbBool)config->enable_hme_level0_flag;
    callback_data->eb_enc_parameters.enable_hme_level1_flag =
        (EbBool)config->enable_hme_level1_flag;
    callback_data->eb_enc_parameters.enable_hme_level2_flag =
        (EbBool)config->enable_hme_level2_flag;
    callback_data->eb_enc_parameters.search_area_width  = config->search_area_width;
    callback_data->eb_enc_parameters.search_area_height = config->search_area_height;
    callback_data->eb_enc_parameters.number_hme_search_region_in_width =
        config->number_hme_search_region_in_width;
    callback_data->eb_enc_parameters.number_hme_search_region_in_height =
        config->number_hme_search_region_in_height;
    callback_data->eb_enc_parameters.hme_level0_total_search_area_width =
        config->hme_level0_total_search_area_width;
    callback_data->eb_enc_parameters.hme_level0_total_search_area_height =
        config->hme_level0_total_search_area_height;
    callback_data->eb_enc_parameters.screen_content_mode = (EbBool)config->screen_content_mode;
    callback_data->eb_enc_parameters.intrabc_mode = config->intrabc_mode;
    callback_data->eb_enc_parameters.enable_hbd_mode_decision =
        config->enable_hbd_mode_decision;
    callback_data->eb_enc_parameters.palette_level            = config->palette_level;
    callback_data->eb_enc_parameters.channel_id               = config->channel_id;
    callback_data->eb_enc_parameters.active_channel_count     = config->active_channel_count;
    callback_data->eb_enc_parameters.high_dynamic_range_input = config->high_dynamic_range_input;
    callback_data->eb_enc_parameters.encoder_bit_depth        = config->encoder_bit_depth;
    callback_data->eb_enc_parameters.is_16bit_pipeline   = config->is_16bit_pipeline;
    callback_data->eb_enc_parameters.encoder_color_format =
        (EbColorFormat)config->encoder_color_format;
    callback_data->eb_enc_parameters.compressed_ten_bit_format = config->compressed_ten_bit_format;
    callback_data->eb_enc_parameters.profile                   = config->profile;
    callback_data->eb_enc_parameters.tier                      = config->tier;
    callback_data->eb_enc_parameters.level                     = config->level;
    callback_data->eb_enc_parameters.injector_frame_rate       = config->injector_frame_rate;
    callback_data->eb_enc_parameters.speed_control_flag        = config->speed_control_flag;
    callback_data->eb_enc_parameters.use_cpu_flags             = config->cpu_flags_limit;
    callback_data->eb_enc_parameters.logical_processors        = config->logical_processors;
    callback_data->eb_enc_parameters.unpin                 = config->unpin;
    callback_data->eb_enc_parameters.target_socket             = config->target_socket;
    callback_data->eb_enc_parameters.unrestricted_motion_vector =
        config->unrestricted_motion_vector;
    callback_data->eb_enc_parameters.recon_enabled = config->recon_file ? EB_TRUE : EB_FALSE;
    // --- start: ALTREF_FILTERING_SUPPORT
    callback_data->eb_enc_parameters.tf_level = (int8_t)config->tf_level;
    callback_data->eb_enc_parameters.altref_strength = config->altref_strength;
    callback_data->eb_enc_parameters.altref_nframes  = config->altref_nframes;
    callback_data->eb_enc_parameters.enable_overlays = (EbBool)config->enable_overlays;
    // --- end: ALTREF_FILTERING_SUPPORT
    // --- start: SUPER-RESOLUTION SUPPORT
    callback_data->eb_enc_parameters.superres_mode = config->superres_mode;
    callback_data->eb_enc_parameters.superres_denom = config->superres_denom;
    callback_data->eb_enc_parameters.superres_kf_denom = config->superres_kf_denom;
    callback_data->eb_enc_parameters.superres_qthres = config->superres_qthres;
    // --- end: SUPER-RESOLUTION SUPPORT

    for (hme_region_index = 0;
         hme_region_index < callback_data->eb_enc_parameters.number_hme_search_region_in_width;
         ++hme_region_index) {
        callback_data->eb_enc_parameters.hme_level0_search_area_in_width_array[hme_region_index] =
            config->hme_level0_search_area_in_width_array[hme_region_index];
        callback_data->eb_enc_parameters.hme_level1_search_area_in_width_array[hme_region_index] =
            config->hme_level1_search_area_in_width_array[hme_region_index];
        callback_data->eb_enc_parameters.hme_level2_search_area_in_width_array[hme_region_index] =
            config->hme_level2_search_area_in_width_array[hme_region_index];
    }

    for (hme_region_index = 0;
         hme_region_index < callback_data->eb_enc_parameters.number_hme_search_region_in_height;
         ++hme_region_index) {
        callback_data->eb_enc_parameters.hme_level0_search_area_in_height_array[hme_region_index] =
            config->hme_level0_search_area_in_height_array[hme_region_index];
        callback_data->eb_enc_parameters.hme_level1_search_area_in_height_array[hme_region_index] =
            config->hme_level1_search_area_in_height_array[hme_region_index];
        callback_data->eb_enc_parameters.hme_level2_search_area_in_height_array[hme_region_index] =
            config->hme_level2_search_area_in_height_array[hme_region_index];
    }

    // Prediction Structure
    callback_data->eb_enc_parameters.enable_manual_pred_struct    = config->enable_manual_pred_struct;
    if (callback_data->eb_enc_parameters.enable_manual_pred_struct) {
        callback_data->eb_enc_parameters.manual_pred_struct_entry_num = config->manual_pred_struct_entry_num;
        memcpy(&callback_data->eb_enc_parameters.pred_struct[0], &config->pred_struct[0],config->manual_pred_struct_entry_num*sizeof(PredictionStructureConfigEntry));
    }

    return return_error;
}

static EbErrorType allocate_frame_buffer(EbConfig *config, uint8_t *p_buffer) {
    const int32_t ten_bit_packed_mode =
        (config->encoder_bit_depth > 8) && (config->compressed_ten_bit_format == 0) ? 1 : 0;

    // Chroma subsampling
    const EbColorFormat color_format  = (EbColorFormat)config->encoder_color_format;
    const uint8_t       subsampling_x = (color_format == EB_YUV444 ? 1 : 2) - 1;

    // Determine size of each plane
    const size_t luma_8bit_size =
        config->input_padded_width * config->input_padded_height * (1 << ten_bit_packed_mode);

    const size_t chroma_8bit_size = luma_8bit_size >> (3 - color_format);

    const size_t luma_10bit_size =
        (config->encoder_bit_depth > 8 && ten_bit_packed_mode == 0) ? luma_8bit_size : 0;

    const size_t chroma_10bit_size =
        (config->encoder_bit_depth > 8 && ten_bit_packed_mode == 0) ? chroma_8bit_size : 0;

    // Determine
    EbSvtIOFormat *input_ptr = (EbSvtIOFormat *)p_buffer;
    input_ptr->y_stride      = config->input_padded_width;
    input_ptr->cr_stride     = config->input_padded_width >> subsampling_x;
    input_ptr->cb_stride     = config->input_padded_width >> subsampling_x;

    if (luma_8bit_size) {
        EB_APP_MALLOC(
            uint8_t *, input_ptr->luma, luma_8bit_size, EB_N_PTR, EB_ErrorInsufficientResources);
    } else {
        input_ptr->luma = 0;
    }

    if (chroma_8bit_size) {
        EB_APP_MALLOC(
            uint8_t *, input_ptr->cb, chroma_8bit_size, EB_N_PTR, EB_ErrorInsufficientResources);
    } else {
        input_ptr->cb = 0;
    }

    if (chroma_8bit_size) {
        EB_APP_MALLOC(
            uint8_t *, input_ptr->cr, chroma_8bit_size, EB_N_PTR, EB_ErrorInsufficientResources);
    } else {
        input_ptr->cr = 0;
    }

    if (luma_10bit_size) {
        EB_APP_MALLOC(uint8_t *,
                      input_ptr->luma_ext,
                      luma_10bit_size,
                      EB_N_PTR,
                      EB_ErrorInsufficientResources);
    } else {
        input_ptr->luma_ext = 0;
    }

    if (chroma_10bit_size) {
        EB_APP_MALLOC(uint8_t *,
                      input_ptr->cb_ext,
                      chroma_10bit_size,
                      EB_N_PTR,
                      EB_ErrorInsufficientResources);
    } else {
        input_ptr->cb_ext = 0;
    }

    if (chroma_10bit_size) {
        EB_APP_MALLOC(uint8_t *,
                      input_ptr->cr_ext,
                      chroma_10bit_size,
                      EB_N_PTR,
                      EB_ErrorInsufficientResources);
    } else {
        input_ptr->cr_ext = 0;
    }

    return EB_ErrorNone;
}

EbErrorType allocate_input_buffers(EbConfig *config, EbAppContext *callback_data) {
    EB_APP_MALLOC(EbBufferHeaderType *,
                  callback_data->input_buffer_pool,
                  sizeof(EbBufferHeaderType),
                  EB_N_PTR,
                  EB_ErrorInsufficientResources);

    // Initialize Header
    callback_data->input_buffer_pool->size = sizeof(EbBufferHeaderType);

    EB_APP_MALLOC(uint8_t *,
                  callback_data->input_buffer_pool->p_buffer,
                  sizeof(EbSvtIOFormat),
                  EB_N_PTR,
                  EB_ErrorInsufficientResources);

    // Allocate frame buffer for the p_buffer
    if (config->buffered_input == -1)
        allocate_frame_buffer(config, callback_data->input_buffer_pool->p_buffer);

    // Assign the variables
    callback_data->input_buffer_pool->p_app_private = NULL;
    callback_data->input_buffer_pool->pic_type      = EB_AV1_INVALID_PICTURE;

    return EB_ErrorNone;
}

EbErrorType allocate_output_recon_buffers(EbConfig *config, EbAppContext *callback_data) {
    const size_t luma_size = config->input_padded_width * config->input_padded_height;

    // both u and v
    const size_t chroma_size = luma_size >> (3 - config->encoder_color_format);
    const size_t ten_bit     = (config->encoder_bit_depth > 8);
    const size_t frame_size  = (luma_size + 2 * chroma_size) << ten_bit;

    // Recon Port
    EB_APP_MALLOC(EbBufferHeaderType *,
                  callback_data->recon_buffer,
                  sizeof(EbBufferHeaderType),
                  EB_N_PTR,
                  EB_ErrorInsufficientResources);

    // Initialize Header
    callback_data->recon_buffer->size = sizeof(EbBufferHeaderType);

    EB_APP_MALLOC(uint8_t *,
                  callback_data->recon_buffer->p_buffer,
                  frame_size,
                  EB_N_PTR,
                  EB_ErrorInsufficientResources);

    callback_data->recon_buffer->n_alloc_len   = (uint32_t)frame_size;
    callback_data->recon_buffer->p_app_private = NULL;

    return EB_ErrorNone;
}

EbErrorType preload_frames_info_ram(EbConfig *config) {
    EbErrorType         return_error = EB_ErrorNone;
    int32_t             input_padded_width  = config->input_padded_width;
    int32_t             input_padded_height = config->input_padded_height;
    size_t              read_size;
    const EbColorFormat color_format = (EbColorFormat)
                                           config->encoder_color_format; // Chroma subsampling

    read_size = input_padded_width * input_padded_height; //Luma
    read_size += 2 * (read_size >> (3 - color_format)); // Add Chroma
    if (config->encoder_bit_depth == 10 && config->compressed_ten_bit_format == 1)
        read_size += read_size / 4;
    else if (config->encoder_bit_depth > 8)
        read_size *= 2; //10 bit
    EB_APP_MALLOC(uint8_t **,
                  config->sequence_buffer,
                  sizeof(uint8_t *) * config->buffered_input,
                  EB_N_PTR,
                  EB_ErrorInsufficientResources);

    for (int32_t processed_frame_count = 0; processed_frame_count < config->buffered_input;
         ++processed_frame_count) {
        EB_APP_MALLOC(uint8_t *,
                      config->sequence_buffer[processed_frame_count],
                      read_size,
                      EB_N_PTR,
                      EB_ErrorInsufficientResources);
        // Fill the buffer with a complete frame
        size_t filled_len = fread(
            config->sequence_buffer[processed_frame_count], 1, read_size, config->input_file);

        if (read_size != filled_len) {
            fseek(config->input_file, 0, SEEK_SET);

            // Fill the buffer with a complete frame
            if (read_size != fread(config->sequence_buffer[processed_frame_count], 1, read_size, config->input_file))
                return_error = EB_Corrupt_Frame;
        }
    }

    return return_error;
}

/***************************************
* Functions Implementation
***************************************/

/***********************************
 * Initialize Core & Component
 ***********************************/
EbErrorType init_encoder(EbConfig *config, EbAppContext *callback_data, uint32_t instance_idx) {
    // Allocate a memory table hosting all allocated pointers
    allocate_memory_table(instance_idx);

    ///************************* LIBRARY INIT [START] *********************///
    // STEP 1: Call the library to construct a Component Handle
    EbErrorType return_error = svt_av1_enc_init_handle(
        &callback_data->svt_encoder_handle, callback_data, &callback_data->eb_enc_parameters);

    if (return_error != EB_ErrorNone) return return_error;
    // STEP 3: Copy all configuration parameters into the callback structure
    return_error = copy_configuration_parameters(config, callback_data, instance_idx);

    if (return_error != EB_ErrorNone) return return_error;
    // STEP 4: Send over all configuration parameters
    // Set the Parameters
    return_error = svt_av1_enc_set_parameter(callback_data->svt_encoder_handle,
                                            &callback_data->eb_enc_parameters);

    if (return_error != EB_ErrorNone) return return_error;
    // STEP 5: Init Encoder
    return_error = svt_av1_enc_init(callback_data->svt_encoder_handle);
    if (return_error != EB_ErrorNone) { return return_error; }

    ///************************* LIBRARY INIT [END] *********************///

    ///********************** APPLICATION INIT [START] ******************///

    // STEP 6: Allocate input buffers carrying the yuv frames in
    return_error = allocate_input_buffers(config, callback_data);

    if (return_error != EB_ErrorNone) return return_error;
    // STEP 7: Allocate output Recon Buffer
    return_error = allocate_output_recon_buffers(config, callback_data);

    if (return_error != EB_ErrorNone) return return_error;
    // Allocate the Sequence Buffer
    if (config->buffered_input != -1) {
        // Preload frames into the ram for a faster yuv access time
        preload_frames_info_ram(config);
    } else
        config->sequence_buffer = 0;
    ///********************** APPLICATION INIT [END] ******************////////

    return return_error;
}

/***********************************
 * Deinit Components
 ***********************************/
EbErrorType de_init_encoder(EbAppContext *callback_data_ptr, uint32_t instance_index) {
    EbErrorType       return_error = EB_ErrorNone;
    int32_t           ptr_index    = 0;
    EbMemoryMapEntry *memory_entry = (EbMemoryMapEntry *)0;

    if (((EbComponentType *)(callback_data_ptr->svt_encoder_handle)) != NULL)
        return_error = svt_av1_enc_deinit(callback_data_ptr->svt_encoder_handle);
    // Destruct the buffer memory pool
    if (return_error != EB_ErrorNone) return return_error;
    // Loop through the ptr table and free all malloc'd pointers per channel
    for (ptr_index = app_memory_map_index_all_channels[instance_index] - 1; ptr_index >= 0;
         --ptr_index) {
        memory_entry = &app_memory_map_all_channels[instance_index][ptr_index];
        switch (memory_entry->ptr_type) {
        case EB_N_PTR: free(memory_entry->ptr); break;
        default: return_error = EB_ErrorMax; break;
        }
    }
    free(app_memory_map_all_channels[instance_index]);

    // Destruct the component
    svt_av1_enc_deinit_handle(callback_data_ptr->svt_encoder_handle);

    return return_error;
}
