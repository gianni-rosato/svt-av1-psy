/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

/***************************************
 * Includes
 ***************************************/

#include <stdlib.h>

#include "EbAppContext.h"
#include "EbAppConfig.h"

#define INPUT_SIZE_576p_TH                0x90000        // 0.58 Million
#define INPUT_SIZE_1080i_TH                0xB71B0        // 0.75 Million
#define INPUT_SIZE_1080p_TH                0x1AB3F0    // 1.75 Million
#define INPUT_SIZE_4K_TH                0x29F630    // 2.75 Million

#define IS_16_BIT(bit_depth) (bit_depth==10?1:0)
#define EB_OUTPUTSTREAMBUFFERSIZE_MACRO(ResolutionSize)                ((ResolutionSize) < (INPUT_SIZE_1080i_TH) ? 0x1E8480 : (ResolutionSize) < (INPUT_SIZE_1080p_TH) ? 0x2DC6C0 : (ResolutionSize) < (INPUT_SIZE_4K_TH) ? 0x2DC6C0 : 0x2DC6C0  )

 /***************************************
 * Variables Defining a memory table
 *  hosting all allocated pointers
 ***************************************/
EbMemoryMapEntry                 *app_memory_map;
uint32_t                         *app_memory_map_index;
uint64_t                         *total_app_memory;
uint32_t                          app_malloc_count = 0;
static EbMemoryMapEntry          *appMemoryMapAllChannels[MAX_CHANNEL_NUMBER];
static uint32_t                   appMemoryMapIndexAllChannels[MAX_CHANNEL_NUMBER];
static uint64_t                   appMemoryMallocdAllChannels[MAX_CHANNEL_NUMBER];

/***************************************
* Allocation and initializing a memory table
*  hosting all allocated pointers
***************************************/
void AllocateMemoryTable(
    uint32_t    instance_idx)
{
    // Malloc Memory Table for the instance @ instance_idx
    appMemoryMapAllChannels[instance_idx]        = (EbMemoryMapEntry*)malloc(sizeof(EbMemoryMapEntry) * MAX_APP_NUM_PTR);

    // Init the table index
    appMemoryMapIndexAllChannels[instance_idx]   = 0;

    // Size of the table
    appMemoryMallocdAllChannels[instance_idx]    = sizeof(EbMemoryMapEntry) * MAX_APP_NUM_PTR;
    total_app_memory = &appMemoryMallocdAllChannels[instance_idx];

    // Set pointer to the first entry
    app_memory_map                                = appMemoryMapAllChannels[instance_idx];

    // Set index to the first entry
    app_memory_map_index                           = &appMemoryMapIndexAllChannels[instance_idx];

    // Init Number of pointers
    app_malloc_count = 0;

    return;
}

/*************************************
**************************************
*** Helper functions Input / Output **
**************************************
**************************************/
/******************************************************
* Copy fields from the stream to the input buffer
    Input   : stream
    Output  : valid input buffer
******************************************************/
void ProcessInputFieldBufferingMode(
    uint64_t                      processed_frame_count,
    int32_t                      *filledLen,
    FILE                         *input_file,
    uint8_t                        *lumaInputPtr,
    uint8_t                        *cbInputPtr,
    uint8_t                        *crInputPtr,
    uint32_t                      input_padded_width,
    uint32_t                      input_padded_height,
    uint8_t                       is16bit) {
    uint64_t  sourceLumaRowSize   = (uint64_t)(input_padded_width << is16bit);
    uint64_t  sourceChromaRowSize = sourceLumaRowSize >> 1;

    uint8_t  *ebInputPtr;
    uint32_t  inputRowIndex;

    // Y
    ebInputPtr = lumaInputPtr;
    // Skip 1 luma row if bottom field (point to the bottom field)
    if (processed_frame_count % 2 != 0)
        fseeko(input_file, (long)sourceLumaRowSize, SEEK_CUR);

    for (inputRowIndex = 0; inputRowIndex < input_padded_height; inputRowIndex++) {
        *filledLen += (uint32_t)fread(ebInputPtr, 1, sourceLumaRowSize, input_file);
        // Skip 1 luma row (only fields)
        fseeko(input_file, (long)sourceLumaRowSize, SEEK_CUR);
        ebInputPtr += sourceLumaRowSize;
    }

    // U
    ebInputPtr = cbInputPtr;
    // Step back 1 luma row if bottom field (undo the previous jump), and skip 1 chroma row if bottom field (point to the bottom field)
    if (processed_frame_count % 2 != 0) {
        fseeko(input_file, -(long)sourceLumaRowSize, SEEK_CUR);
        fseeko(input_file, (long)sourceChromaRowSize, SEEK_CUR);
    }

    for (inputRowIndex = 0; inputRowIndex < input_padded_height >> 1; inputRowIndex++) {
        *filledLen += (uint32_t)fread(ebInputPtr, 1, sourceChromaRowSize, input_file);
        // Skip 1 chroma row (only fields)
        fseeko(input_file, (long)sourceChromaRowSize, SEEK_CUR);
        ebInputPtr += sourceChromaRowSize;
    }

    // V
    ebInputPtr = crInputPtr;
    // Step back 1 chroma row if bottom field (undo the previous jump), and skip 1 chroma row if bottom field (point to the bottom field)
    // => no action

    for (inputRowIndex = 0; inputRowIndex < input_padded_height >> 1; inputRowIndex++) {
        *filledLen += (uint32_t)fread(ebInputPtr, 1, sourceChromaRowSize, input_file);
        // Skip 1 chroma row (only fields)
        fseeko(input_file, (long)sourceChromaRowSize, SEEK_CUR);
        ebInputPtr += sourceChromaRowSize;
    }

    // Step back 1 chroma row if bottom field (undo the previous jump)
    if (processed_frame_count % 2 != 0)
        fseeko(input_file, -(long)sourceChromaRowSize, SEEK_CUR);
}

/***********************************************
* Copy configuration parameters from
*  The config structure, to the
*  callback structure to send to the library
***********************************************/
EbErrorType CopyConfigurationParameters(
    EbConfig                *config,
    EbAppContext            *callback_data,
    uint32_t                 instance_idx)
{
    EbErrorType   return_error = EB_ErrorNone;
    uint32_t         hmeRegionIndex;

    // Assign Instance index to the library
    callback_data->instance_idx = (uint8_t)instance_idx;

    // Initialize Port Activity Flags
    callback_data->output_stream_port_active = APP_PortActive;

    callback_data->eb_enc_parameters.source_width = config->source_width;
    callback_data->eb_enc_parameters.source_height = config->source_height;
    callback_data->eb_enc_parameters.render_width = config->input_padded_width;
    callback_data->eb_enc_parameters.render_height = config->input_padded_height;

    callback_data->eb_enc_parameters.intra_period_length = config->intra_period;
    callback_data->eb_enc_parameters.intra_refresh_type = config->intra_refresh_type;
    callback_data->eb_enc_parameters.base_layer_switch_mode = config->base_layer_switch_mode;
    callback_data->eb_enc_parameters.enc_mode = (EbBool)config->enc_mode;
#if TWO_PASS_USE_2NDP_ME_IN_1STP
    callback_data->eb_enc_parameters.snd_pass_enc_mode = (EbBool)config->snd_pass_enc_mode;
#endif
    callback_data->eb_enc_parameters.frame_rate = config->frame_rate;
    callback_data->eb_enc_parameters.frame_rate_denominator = config->frame_rate_denominator;
    callback_data->eb_enc_parameters.frame_rate_numerator = config->frame_rate_numerator;
    callback_data->eb_enc_parameters.hierarchical_levels = config->hierarchical_levels;
    callback_data->eb_enc_parameters.pred_structure = (uint8_t)config->pred_structure;
    callback_data->eb_enc_parameters.in_loop_me_flag = config->in_loop_me_flag;
    callback_data->eb_enc_parameters.ext_block_flag = config->ext_block_flag;
    callback_data->eb_enc_parameters.tile_rows = config->tile_rows;
    callback_data->eb_enc_parameters.tile_columns = config->tile_columns;

    callback_data->eb_enc_parameters.scene_change_detection = config->scene_change_detection;
    callback_data->eb_enc_parameters.look_ahead_distance = config->look_ahead_distance;
    callback_data->eb_enc_parameters.frames_to_be_encoded = config->frames_to_be_encoded;
    callback_data->eb_enc_parameters.rate_control_mode = config->rate_control_mode;
    callback_data->eb_enc_parameters.target_bit_rate = config->target_bit_rate;
    callback_data->eb_enc_parameters.max_qp_allowed = config->max_qp_allowed;
    callback_data->eb_enc_parameters.min_qp_allowed = config->min_qp_allowed;
    callback_data->eb_enc_parameters.enable_adaptive_quantization = (EbBool)config->enable_adaptive_quantization;
    callback_data->eb_enc_parameters.qp = config->qp;
    callback_data->eb_enc_parameters.use_qp_file = (EbBool)config->use_qp_file;
#if TWO_PASS
    callback_data->eb_enc_parameters.input_stat_file = config->input_stat_file;
    callback_data->eb_enc_parameters.output_stat_file = config->output_stat_file;
#endif
    callback_data->eb_enc_parameters.stat_report = (EbBool)config->stat_report;
    callback_data->eb_enc_parameters.disable_dlf_flag = (EbBool)config->disable_dlf_flag;
    callback_data->eb_enc_parameters.enable_warped_motion = (EbBool)config->enable_warped_motion;
    callback_data->eb_enc_parameters.enable_global_motion = (EbBool)config->enable_global_motion;
    callback_data->eb_enc_parameters.enable_obmc = (EbBool)config->enable_obmc;
    callback_data->eb_enc_parameters.enable_rdoq = config->enable_rdoq;
    callback_data->eb_enc_parameters.enable_filter_intra = (EbBool)config->enable_filter_intra;
    callback_data->eb_enc_parameters.use_default_me_hme = (EbBool)config->use_default_me_hme;
    callback_data->eb_enc_parameters.enable_hme_flag = (EbBool)config->enable_hme_flag;
    callback_data->eb_enc_parameters.enable_hme_level0_flag = (EbBool)config->enable_hme_level0_flag;
    callback_data->eb_enc_parameters.enable_hme_level1_flag = (EbBool)config->enable_hme_level1_flag;
    callback_data->eb_enc_parameters.enable_hme_level2_flag = (EbBool)config->enable_hme_level2_flag;
    callback_data->eb_enc_parameters.search_area_width = config->search_area_width;
    callback_data->eb_enc_parameters.search_area_height = config->search_area_height;
    callback_data->eb_enc_parameters.number_hme_search_region_in_width = config->number_hme_search_region_in_width;
    callback_data->eb_enc_parameters.number_hme_search_region_in_height = config->number_hme_search_region_in_height;
    callback_data->eb_enc_parameters.hme_level0_total_search_area_width = config->hme_level0_total_search_area_width;
    callback_data->eb_enc_parameters.hme_level0_total_search_area_height = config->hme_level0_total_search_area_height;
    callback_data->eb_enc_parameters.screen_content_mode = (EbBool)config->screen_content_mode;
    callback_data->eb_enc_parameters.enable_hbd_mode_decision = (EbBool)config->enable_hbd_mode_decision;
    callback_data->eb_enc_parameters.enable_palette = config->enable_palette;
    callback_data->eb_enc_parameters.olpd_refinement = config->olpd_refinement;
    callback_data->eb_enc_parameters.constrained_intra = (EbBool)config->constrained_intra;
    callback_data->eb_enc_parameters.channel_id = config->channel_id;
    callback_data->eb_enc_parameters.active_channel_count = config->active_channel_count;
    callback_data->eb_enc_parameters.high_dynamic_range_input = config->high_dynamic_range_input;
    callback_data->eb_enc_parameters.encoder_bit_depth = config->encoder_bit_depth;
    callback_data->eb_enc_parameters.encoder_color_format =
        (EbColorFormat)config->encoder_color_format;
    callback_data->eb_enc_parameters.compressed_ten_bit_format = config->compressed_ten_bit_format;
    callback_data->eb_enc_parameters.profile = config->profile;
    callback_data->eb_enc_parameters.tier = config->tier;
    callback_data->eb_enc_parameters.level = config->level;
    callback_data->eb_enc_parameters.injector_frame_rate = config->injector_frame_rate;
    callback_data->eb_enc_parameters.speed_control_flag = config->speed_control_flag;
    callback_data->eb_enc_parameters.use_cpu_flags = config->cpu_flags_limit;
    callback_data->eb_enc_parameters.logical_processors = config->logical_processors;
    callback_data->eb_enc_parameters.target_socket = config->target_socket;
    callback_data->eb_enc_parameters.unrestricted_motion_vector = config->unrestricted_motion_vector;
    callback_data->eb_enc_parameters.recon_enabled = config->recon_file ? EB_TRUE : EB_FALSE;
    // --- start: ALTREF_FILTERING_SUPPORT
    callback_data->eb_enc_parameters.enable_altrefs  = (EbBool)config->enable_altrefs;
    callback_data->eb_enc_parameters.altref_strength = config->altref_strength;
    callback_data->eb_enc_parameters.altref_nframes  = config->altref_nframes;
    callback_data->eb_enc_parameters.enable_overlays = (EbBool)config->enable_overlays;
    // --- end: ALTREF_FILTERING_SUPPORT

    for (hmeRegionIndex = 0; hmeRegionIndex < callback_data->eb_enc_parameters.number_hme_search_region_in_width; ++hmeRegionIndex) {
        callback_data->eb_enc_parameters.hme_level0_search_area_in_width_array[hmeRegionIndex] = config->hme_level0_search_area_in_width_array[hmeRegionIndex];
        callback_data->eb_enc_parameters.hme_level1_search_area_in_width_array[hmeRegionIndex] = config->hme_level1_search_area_in_width_array[hmeRegionIndex];
        callback_data->eb_enc_parameters.hme_level2_search_area_in_width_array[hmeRegionIndex] = config->hme_level2_search_area_in_width_array[hmeRegionIndex];
    }

    for (hmeRegionIndex = 0; hmeRegionIndex < callback_data->eb_enc_parameters.number_hme_search_region_in_height; ++hmeRegionIndex) {
        callback_data->eb_enc_parameters.hme_level0_search_area_in_height_array[hmeRegionIndex] = config->hme_level0_search_area_in_height_array[hmeRegionIndex];
        callback_data->eb_enc_parameters.hme_level1_search_area_in_height_array[hmeRegionIndex] = config->hme_level1_search_area_in_height_array[hmeRegionIndex];
        callback_data->eb_enc_parameters.hme_level2_search_area_in_height_array[hmeRegionIndex] = config->hme_level2_search_area_in_height_array[hmeRegionIndex];
    }

    callback_data->eb_enc_parameters.sq_weight = config->sq_weight;
    callback_data->eb_enc_parameters.enable_auto_max_partition = config->enable_auto_max_partition;

    callback_data->eb_enc_parameters.md_stage_1_cand_prune_th = config->md_stage_1_cand_prune_th;
    callback_data->eb_enc_parameters.md_stage_1_class_prune_th = config->md_stage_1_class_prune_th;
    callback_data->eb_enc_parameters.md_stage_2_cand_prune_th = config->md_stage_2_cand_prune_th;
    callback_data->eb_enc_parameters.md_stage_2_class_prune_th = config->md_stage_2_class_prune_th;

    return return_error;
}

static EbErrorType AllocateFrameBuffer(EbConfig *config, uint8_t *p_buffer)
{
    const int32_t tenBitPackedMode =
        (config->encoder_bit_depth > 8) &&
        (config->compressed_ten_bit_format == 0) ? 1 : 0;

    // Chroma subsampling
    const EbColorFormat color_format =
        (EbColorFormat)config->encoder_color_format;
    const uint8_t subsampling_x = (color_format == EB_YUV444 ? 1 : 2) - 1;

    // Determine size of each plane
    const size_t luma8bitSize = config->input_padded_width *
        config->input_padded_height * (1 << tenBitPackedMode);

    const size_t chroma8bitSize = luma8bitSize >> (3 - color_format);

    const size_t luma10bitSize =
        (config->encoder_bit_depth > 8 && tenBitPackedMode == 0) ?
            luma8bitSize : 0;

    const size_t chroma10bitSize =
        (config->encoder_bit_depth > 8 && tenBitPackedMode == 0) ?
            chroma8bitSize : 0;

    // Determine
    EbSvtIOFormat* inputPtr = (EbSvtIOFormat*)p_buffer;
    inputPtr->y_stride = config->input_padded_width;
    inputPtr->cr_stride = config->input_padded_width >> subsampling_x;
    inputPtr->cb_stride = config->input_padded_width >> subsampling_x;

    if (luma8bitSize) {
        EB_APP_MALLOC(uint8_t*, inputPtr->luma, luma8bitSize, EB_N_PTR,
                      EB_ErrorInsufficientResources);
    } else {
        inputPtr->luma = 0;
    }

    if (chroma8bitSize) {
        EB_APP_MALLOC(uint8_t*, inputPtr->cb, chroma8bitSize, EB_N_PTR,
                      EB_ErrorInsufficientResources);
    } else {
        inputPtr->cb = 0;
    }

    if (chroma8bitSize) {
        EB_APP_MALLOC(uint8_t*, inputPtr->cr, chroma8bitSize, EB_N_PTR,
                      EB_ErrorInsufficientResources);
    } else {
        inputPtr->cr = 0;
    }

    if (luma10bitSize) {
        EB_APP_MALLOC(uint8_t*, inputPtr->luma_ext, luma10bitSize, EB_N_PTR,
                      EB_ErrorInsufficientResources);
    } else {
        inputPtr->luma_ext = 0;
    }

    if (chroma10bitSize) {
        EB_APP_MALLOC(uint8_t*, inputPtr->cb_ext, chroma10bitSize, EB_N_PTR,
                      EB_ErrorInsufficientResources);
    } else {
        inputPtr->cb_ext = 0;
    }

    if (chroma10bitSize) {
        EB_APP_MALLOC(uint8_t*, inputPtr->cr_ext, chroma10bitSize, EB_N_PTR,
                      EB_ErrorInsufficientResources);
    } else {
        inputPtr->cr_ext = 0;
    }

    return EB_ErrorNone;
}

EbErrorType AllocateInputBuffers(EbConfig *config, EbAppContext *callback_data)
{
    EB_APP_MALLOC(EbBufferHeaderType*, callback_data->input_buffer_pool,
                  sizeof(EbBufferHeaderType), EB_N_PTR,
                  EB_ErrorInsufficientResources);

    // Initialize Header
    callback_data->input_buffer_pool->size = sizeof(EbBufferHeaderType);

    EB_APP_MALLOC(uint8_t*, callback_data->input_buffer_pool->p_buffer,
                  sizeof(EbSvtIOFormat), EB_N_PTR,
                  EB_ErrorInsufficientResources);

    // Allocate frame buffer for the p_buffer
    if (config->buffered_input == -1)
        AllocateFrameBuffer(config, callback_data->input_buffer_pool->p_buffer);

    // Assign the variables
    callback_data->input_buffer_pool->p_app_private = NULL;
    callback_data->input_buffer_pool->pic_type = EB_AV1_INVALID_PICTURE;

    return EB_ErrorNone;
}

EbErrorType AllocateOutputReconBuffers(EbConfig *config,
                                       EbAppContext *callback_data)
{
    const size_t luma_size =
        config->input_padded_width * config->input_padded_height;

    // both u and v
    const size_t chromaSize = luma_size >> (3 - config->encoder_color_format);
    const size_t tenBit = (config->encoder_bit_depth > 8);
    const size_t frameSize = (luma_size + 2 * chromaSize) << tenBit;

    // Recon Port
    EB_APP_MALLOC(EbBufferHeaderType*, callback_data->recon_buffer,
                  sizeof(EbBufferHeaderType), EB_N_PTR,
                  EB_ErrorInsufficientResources);

    // Initialize Header
    callback_data->recon_buffer->size = sizeof(EbBufferHeaderType);

    EB_APP_MALLOC(uint8_t*, callback_data->recon_buffer->p_buffer, frameSize,
                  EB_N_PTR, EB_ErrorInsufficientResources);

    callback_data->recon_buffer->n_alloc_len = (uint32_t)frameSize;
    callback_data->recon_buffer->p_app_private = NULL;

    return EB_ErrorNone;
}

EbErrorType AllocateOutputBuffers(EbConfig *config, EbAppContext *callback_data)
{
    uint32_t outputStreamBufferSize =
        (uint32_t)(EB_OUTPUTSTREAMBUFFERSIZE_MACRO(
                    config->input_padded_height * config->input_padded_width));

    EB_APP_MALLOC(EbBufferHeaderType*, callback_data->stream_buffer_pool,
                  sizeof(EbBufferHeaderType), EB_N_PTR,
                  EB_ErrorInsufficientResources);

    // Initialize Header
    callback_data->stream_buffer_pool->size = sizeof(EbBufferHeaderType);

    EB_APP_MALLOC(uint8_t*, callback_data->stream_buffer_pool->p_buffer,
                  outputStreamBufferSize, EB_N_PTR,
                  EB_ErrorInsufficientResources);

    callback_data->stream_buffer_pool->n_alloc_len = outputStreamBufferSize;
    callback_data->stream_buffer_pool->p_app_private = NULL;
    callback_data->stream_buffer_pool->pic_type = EB_AV1_INVALID_PICTURE;

    return EB_ErrorNone;
}

EbErrorType PreloadFramesIntoRam(
    EbConfig                *config)
{
    EbErrorType    return_error = EB_ErrorNone;
    int32_t             processed_frame_count;
    int32_t             filledLen;
    int32_t             input_padded_width = config->input_padded_width;
    int32_t             input_padded_height = config->input_padded_height;
    int32_t             readSize;
    const EbColorFormat color_format =
        (EbColorFormat)config->encoder_color_format;  // Chroma subsampling

    FILE *input_file = config->input_file;

    readSize = input_padded_width * input_padded_height; //Luma
    readSize += 2 * (readSize >> (3 - color_format)); // Add Chroma
    if (config->encoder_bit_depth == 10 && config->compressed_ten_bit_format == 1) {
        readSize += readSize / 4;
    } else
        readSize *= (config->encoder_bit_depth > 8 ? 2 : 1); //10 bit
    EB_APP_MALLOC(uint8_t **, config->sequence_buffer, sizeof(uint8_t*) * config->buffered_input, EB_N_PTR, EB_ErrorInsufficientResources);

    for (processed_frame_count = 0; processed_frame_count < config->buffered_input; ++processed_frame_count) {
        EB_APP_MALLOC(uint8_t*, config->sequence_buffer[processed_frame_count], readSize, EB_N_PTR, EB_ErrorInsufficientResources);
        // Interlaced Video
        if (config->separate_fields) {
            EbBool is16bit = config->encoder_bit_depth > 8;
            if (is16bit == 0 || (is16bit == 1 && config->compressed_ten_bit_format == 0)) {
                const int32_t tenBitPackedMode = (config->encoder_bit_depth > 8) && (config->compressed_ten_bit_format == 0) ? 1 : 0;

                const size_t luma8bitSize =

                    (config->input_padded_width) *
                    (config->input_padded_height) *

                    (1 << tenBitPackedMode);

                const size_t chroma8bitSize = luma8bitSize >> 2;

                filledLen = 0;

                ProcessInputFieldBufferingMode(
                    processed_frame_count,
                    &filledLen,
                    input_file,
                    config->sequence_buffer[processed_frame_count],
                    config->sequence_buffer[processed_frame_count] + luma8bitSize,
                    config->sequence_buffer[processed_frame_count] + luma8bitSize + chroma8bitSize,
                    (uint32_t)input_padded_width,
                    (uint32_t)input_padded_height,
                    is16bit);

                if (readSize != filledLen) {
                    fseek(input_file, 0, SEEK_SET);
                    filledLen = 0;

                    ProcessInputFieldBufferingMode(
                        processed_frame_count,
                        &filledLen,
                        input_file,
                        config->sequence_buffer[processed_frame_count],
                        config->sequence_buffer[processed_frame_count] + luma8bitSize,
                        config->sequence_buffer[processed_frame_count] + luma8bitSize + chroma8bitSize,
                        (uint32_t)input_padded_width,
                        (uint32_t)input_padded_height,
                        is16bit);
                }

                // Reset the pointer position after a top field
                if (processed_frame_count % 2 == 0)
                    fseek(input_file, -(long)(readSize << 1), SEEK_CUR);
            }
            // Unpacked 10 bit
            else {
                const int32_t tenBitPackedMode = (config->encoder_bit_depth > 8) && (config->compressed_ten_bit_format == 0) ? 1 : 0;

                const size_t luma8bitSize =
                    (config->input_padded_width) *
                    (config->input_padded_height) *
                    (1 << tenBitPackedMode);

                const size_t chroma8bitSize = luma8bitSize >> 2;

                const size_t luma10bitSize = (config->encoder_bit_depth > 8 && tenBitPackedMode == 0) ? luma8bitSize : 0;
                const size_t chroma10bitSize = (config->encoder_bit_depth > 8 && tenBitPackedMode == 0) ? chroma8bitSize : 0;

                filledLen = 0;

                ProcessInputFieldBufferingMode(
                    processed_frame_count,
                    &filledLen,
                    input_file,
                    config->sequence_buffer[processed_frame_count],
                    config->sequence_buffer[processed_frame_count] + luma8bitSize,
                    config->sequence_buffer[processed_frame_count] + luma8bitSize + chroma8bitSize,
                    (uint32_t)input_padded_width,
                    (uint32_t)input_padded_height,
                    0);

                ProcessInputFieldBufferingMode(
                    processed_frame_count,
                    &filledLen,
                    input_file,
                    config->sequence_buffer[processed_frame_count] + luma8bitSize + (chroma8bitSize << 1),
                    config->sequence_buffer[processed_frame_count] + luma8bitSize + (chroma8bitSize << 1) + luma10bitSize,
                    config->sequence_buffer[processed_frame_count] + luma8bitSize + (chroma8bitSize << 1) + luma10bitSize + chroma10bitSize,
                    (uint32_t)input_padded_width,
                    (uint32_t)input_padded_height,
                    0);

                if (readSize != filledLen) {
                    fseek(input_file, 0, SEEK_SET);
                    filledLen = 0;

                    ProcessInputFieldBufferingMode(
                        processed_frame_count,
                        &filledLen,
                        input_file,
                        config->sequence_buffer[processed_frame_count],
                        config->sequence_buffer[processed_frame_count] + luma8bitSize,
                        config->sequence_buffer[processed_frame_count] + luma8bitSize + chroma8bitSize,
                        (uint32_t)input_padded_width,
                        (uint32_t)input_padded_height,
                        0);

                    ProcessInputFieldBufferingMode(
                        processed_frame_count,
                        &filledLen,
                        input_file,
                        config->sequence_buffer[processed_frame_count] + luma8bitSize + (chroma8bitSize << 1),
                        config->sequence_buffer[processed_frame_count] + luma8bitSize + (chroma8bitSize << 1) + luma10bitSize,
                        config->sequence_buffer[processed_frame_count] + luma8bitSize + (chroma8bitSize << 1) + luma10bitSize + chroma10bitSize,
                        (uint32_t)input_padded_width,
                        (uint32_t)input_padded_height,
                        0);
                }

                // Reset the pointer position after a top field
                if (processed_frame_count % 2 == 0)
                    fseek(input_file, -(long)(readSize << 1), SEEK_CUR);
            }
        } else {
            // Fill the buffer with a complete frame
            filledLen = 0;
            filledLen += (uint32_t)fread(config->sequence_buffer[processed_frame_count], 1, readSize, input_file);

            if (readSize != filledLen) {
                fseek(config->input_file, 0, SEEK_SET);

                // Fill the buffer with a complete frame
                filledLen = 0;
                filledLen += (uint32_t)fread(config->sequence_buffer[processed_frame_count], 1, readSize, input_file);
            }
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
EbErrorType init_encoder(
    EbConfig              *config,
    EbAppContext          *callback_data,
    uint32_t                 instance_idx)
{
    EbErrorType        return_error = EB_ErrorNone;

    // Allocate a memory table hosting all allocated pointers
    AllocateMemoryTable(instance_idx);

    ///************************* LIBRARY INIT [START] *********************///
    // STEP 1: Call the library to construct a Component Handle
    return_error = eb_init_handle(&callback_data->svt_encoder_handle, callback_data, &callback_data->eb_enc_parameters);

    if (return_error != EB_ErrorNone)
        return return_error;
    // STEP 3: Copy all configuration parameters into the callback structure
    return_error = CopyConfigurationParameters(
                    config,
                    callback_data,
                    instance_idx);

    if (return_error != EB_ErrorNone)
        return return_error;
    // STEP 4: Send over all configuration parameters
    // Set the Parameters
    return_error = eb_svt_enc_set_parameter(
                       callback_data->svt_encoder_handle,
                       &callback_data->eb_enc_parameters);

    if (return_error != EB_ErrorNone)
        return return_error;
    // STEP 5: Init Encoder
    return_error = eb_init_encoder(callback_data->svt_encoder_handle);
    if (return_error != EB_ErrorNone) { return return_error; }

    ///************************* LIBRARY INIT [END] *********************///

    ///********************** APPLICATION INIT [START] ******************///

    // STEP 6: Allocate input buffers carrying the yuv frames in
    return_error = AllocateInputBuffers(config, callback_data);

    if (return_error != EB_ErrorNone)
        return return_error;
    // STEP 7: Allocate output buffers carrying the bitstream out
    return_error = AllocateOutputBuffers(
        config,
        callback_data);

    if (return_error != EB_ErrorNone)
        return return_error;
    // STEP 8: Allocate output Recon Buffer
    return_error = AllocateOutputReconBuffers(
        config,
        callback_data);

    if (return_error != EB_ErrorNone)
        return return_error;
    // Allocate the Sequence Buffer
    if (config->buffered_input != -1) {
        // Preload frames into the ram for a faster yuv access time
        PreloadFramesIntoRam(
            config);
    }
    else
        config->sequence_buffer = 0;
    if (return_error != EB_ErrorNone)
        return return_error;
    ///********************** APPLICATION INIT [END] ******************////////

    return return_error;
}

/***********************************
 * Deinit Components
 ***********************************/
EbErrorType de_init_encoder(
    EbAppContext *callback_data_ptr,
    uint32_t        instance_index)
{
    EbErrorType return_error = EB_ErrorNone;
    int32_t              ptrIndex        = 0;
    EbMemoryMapEntry*   memory_entry     = (EbMemoryMapEntry*)0;

    if (((EbComponentType*)(callback_data_ptr->svt_encoder_handle)) != NULL)
            return_error = eb_deinit_encoder(callback_data_ptr->svt_encoder_handle);
    // Destruct the buffer memory pool
    if (return_error != EB_ErrorNone)
        return return_error;
    // Loop through the ptr table and free all malloc'd pointers per channel
    for (ptrIndex = appMemoryMapIndexAllChannels[instance_index] - 1; ptrIndex >= 0; --ptrIndex) {
        memory_entry = &appMemoryMapAllChannels[instance_index][ptrIndex];
        switch (memory_entry->ptr_type) {
        case EB_N_PTR:
            free(memory_entry->ptr);
            break;
        default:
            return_error = EB_ErrorMax;
            break;
        }
    }
    free(appMemoryMapAllChannels[instance_index]);

    // Destruct the component
    eb_deinit_handle(callback_data_ptr->svt_encoder_handle);

    return return_error;
}
