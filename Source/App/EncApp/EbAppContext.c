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

static EbErrorType allocate_frame_buffer(EbConfig *config, uint8_t *p_buffer) {
    EbSvtAv1EncConfiguration* cfg = &config->config;
    const int32_t ten_bit_packed_mode =
        (cfg->encoder_bit_depth > 8) && (cfg->compressed_ten_bit_format == 0) ? 1 : 0;

    // Chroma subsampling
    const EbColorFormat color_format  = (EbColorFormat)cfg->encoder_color_format;
    const uint8_t       subsampling_x = (color_format == EB_YUV444 ? 1 : 2) - 1;

    // Determine size of each plane
    const size_t luma_8bit_size =
        config->input_padded_width * config->input_padded_height * (1 << ten_bit_packed_mode);

    const size_t chroma_8bit_size = luma_8bit_size >> (3 - color_format);

    const size_t luma_10bit_size =
        (cfg->encoder_bit_depth > 8 && ten_bit_packed_mode == 0) ? luma_8bit_size : 0;

    const size_t chroma_10bit_size =
        (cfg->encoder_bit_depth > 8 && ten_bit_packed_mode == 0) ? chroma_8bit_size : 0;

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
    const size_t chroma_size = luma_size >> (3 - config->config.encoder_color_format);
    const size_t ten_bit     = (config->config.encoder_bit_depth > 8);
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
                                           config->config.encoder_color_format; // Chroma subsampling

    read_size = input_padded_width * input_padded_height; //Luma
    read_size += 2 * (read_size >> (3 - color_format)); // Add Chroma
    if (config->config.encoder_bit_depth == 10 && config->config.compressed_ten_bit_format == 1)
        read_size += read_size / 4;
    else if (config->config.encoder_bit_depth > 8)
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


    callback_data->instance_idx = (uint8_t)instance_idx;
    // Initialize Port Activity Flags
    callback_data->output_stream_port_active                = APP_PortActive;


    // Send over all configuration parameters
    // Set the Parameters
    EbErrorType return_error = svt_av1_enc_set_parameter(callback_data->svt_encoder_handle,
                                            &config->config);

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
