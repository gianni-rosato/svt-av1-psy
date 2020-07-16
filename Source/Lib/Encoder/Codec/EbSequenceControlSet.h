/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

/*
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at www.aomedia.org/license/software. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at www.aomedia.org/license/patent.
*/

#ifndef EbSequenceControlSet_h
#define EbSequenceControlSet_h

#include "EbDefinitions.h"
#include "EbAv1Structs.h"
#include "EbEncodeContext.h"
#include "EbObject.h"

#ifdef __cplusplus
extern "C" {
#endif
/************************************
     * Sequence Control Set
     ************************************/
typedef struct SequenceControlSet {

    /*!< Pointer to the dtor of the struct*/
    EbDctor dctor;
    /*!< Encoding context pointer containing the handle pointer */
    EncodeContext *encode_context_ptr;

    /*!< API structure */
    EbSvtAv1EncConfiguration static_config;
    /*!< Pointer to prediction structure containing the mini-gop information */
    PredictionStructure *    pred_struct_ptr;
    /*!< Super block geomerty pointer */
    SbGeom *                 sb_geom;
    /*!< Array of superblock parameters computed at the resource coordination stage */
    SbParams *               sb_params_array;
    /*!< Bitstream level */
    BitstreamLevel           level[MAX_NUM_OPERATING_POINTS];
    /*!< Sequence header structure, common between the encoder and decoder */
    SeqHeader                seq_header;

    /*!< Sequence coding parameters
            parameters/features are set to be set for the full stream
            but encoding decisions may still be taken at a picture / sub picture level
    */

    /*!< Maximum number of references that a picture can have within the stream needs to be cleaned up*/
    uint32_t max_ref_count;
    /*!< Maximum number of references that a picture can have within the stream */
    uint32_t reference_count;
    /*!< The frequency of intra pictures */
    int32_t  intra_period_length;
    /*!< Intra refresh type 2= key frame, 1= fwd key frame */
    uint32_t intra_refresh_type;
    /*!< Target bitrate in bits per seconds */
    uint32_t target_bitrate;
    /*!< Maximum number of allowed temporal layers */
    uint32_t max_temporal_layers;
    /*!< Overflow bits used for the picture order count increments */
    uint32_t bits_for_picture_order_count;
    /*!< Screen change detection mode 0=OFF, 1=use decimated picture, 2=use full picture */
    EbScdMode scd_mode;
    /*!< Number of delay frames needed to implement future window
         for algorithms such as SceneChange or TemporalFiltering */
    uint32_t scd_delay;
    /*!< Enable the use of altrefs in the stream */
    EbBool    enable_altrefs;
    /*!<  */
    EbBlockMeanPrec          block_mean_calc_prec;
#if !REMOVE_MRP_MODE
    /*!< MRP (The signal changes per preset; 0: MRP mode 0, 1: MRP mode 1) Default is 0. */
    uint8_t mrp_mode;
#endif
    /*!< CDF (The signal changes per preset; 0: CDF update, 1: no CDF update) Default is 0.*/
    uint8_t cdf_mode;
#if !NSQ_REMOVAL_CODE_CLEAN_UP
    /*!< Non-square present flag to use for memory allocation
        (The signal changes per preset; 0: NSQ absent, 1: NSQ present) Default is 1. */
    uint8_t nsq_present;
#endif
    /*!< Down-sampling method @ ME and alt-ref temporal filtering
        (The signal changes per preset; 0: filtering, 1: decimation) Default is 0. */
    uint8_t down_sampling_method_me_search;
    /*!< Allow the usage of motion field motion vector in the stream
        (The signal changes per preset; 0: Enabled, 1: Disabled) Default is 1. */
    uint8_t mfmv_enabled;
    /*!< Film grain strenght */
    int32_t  film_grain_denoise_strength;
    /*!< Film grain seed */
    uint16_t film_grain_random_seed;
    /*!< over_boundary_block: pad resolution to a multiple of SB for smaller overhead
        (The signal changes per preset; 0: No over boundary blk allowed, 1: over boundary blk allowed) Default is 1.
        to enable when md_skip_blk is on */
    uint8_t over_boundary_block_mode;
    /*!< Enable compound prediction to be used in the stream, decisions will be taken at a picture level subsequently
    (The signal changes per preset; 0: compound disabled, 1: compound enabled) Default is 1. */
    uint8_t compound_mode;

    /*!< Temporary input / output statistics files for 2-pass encoding */
    EbBool use_input_stat_file;
    EbBool use_output_stat_file;

    /*!< Sequence resolution parameters */
    uint32_t          chroma_format_idc;
    uint16_t          subsampling_x; // add chroma subsampling parameters
    uint16_t          subsampling_y;
    uint16_t          max_input_luma_width;
    uint16_t          max_input_luma_height;
    uint16_t          max_input_chroma_width;
    uint16_t          max_input_chroma_height;
    uint16_t          max_input_pad_bottom;
    uint16_t          max_input_pad_right;
    uint32_t          chroma_width;
    uint32_t          chroma_height;
    uint32_t          pad_right;
    uint32_t          pad_bottom;
    uint16_t          left_padding;
    uint16_t          top_padding;
    uint16_t          right_padding;
    uint16_t          bot_padding;
    uint32_t          frame_rate;
    uint32_t          encoder_bit_depth;
    EbInputResolution input_resolution;

    /*!< Super block parameters set for the stream */
    uint8_t  sb_sz;
    uint8_t  max_sb_depth;
    uint8_t  pic_width_in_sb;
    uint8_t  picture_height_in_sb;
    uint16_t sb_total_count;
    uint16_t sb_size_pix;
    uint16_t sb_tot_cnt;
    uint16_t max_block_cnt;

    /*!< Block limits */
    uint8_t max_blk_size;
    uint8_t min_blk_size;
    uint8_t max_intra_size;
    uint8_t min_intra_size;

    /*!< Sub picture reagions for picture analysis */
    uint32_t picture_analysis_number_of_regions_per_width;
    uint32_t picture_analysis_number_of_regions_per_height;

    /*!< Tile groups per hierarchical layers */
    uint8_t tile_group_col_count_array[MAX_TEMPORAL_LAYERS];
    uint8_t tile_group_row_count_array[MAX_TEMPORAL_LAYERS];

    /*!< Segements (sub picture) count for different processes */
    uint32_t me_segment_column_count_array[MAX_TEMPORAL_LAYERS];
    uint32_t me_segment_row_count_array[MAX_TEMPORAL_LAYERS];
    uint32_t enc_dec_segment_col_count_array[MAX_TEMPORAL_LAYERS];
    uint32_t enc_dec_segment_row_count_array[MAX_TEMPORAL_LAYERS];
    uint32_t cdef_segment_column_count;
    uint32_t cdef_segment_row_count;
    uint32_t rest_segment_column_count;
    uint32_t rest_segment_row_count;
    uint32_t tf_segment_column_count;
    uint32_t tf_segment_row_count;

    /*!< Picture, reference, recon and input output buffer count */
    uint32_t picture_control_set_pool_init_count;
    uint32_t picture_control_set_pool_init_count_child;
    uint32_t pa_reference_picture_buffer_init_count;
    uint32_t reference_picture_buffer_init_count;
    uint32_t input_buffer_fifo_init_count;
    uint32_t overlay_input_picture_buffer_init_count;
    uint32_t output_stream_buffer_fifo_init_count;
    uint32_t output_recon_buffer_fifo_init_count;

    /*!< Inter processes fifos count */
    uint32_t resource_coordination_fifo_init_count;
    uint32_t picture_analysis_fifo_init_count;
    uint32_t picture_decision_fifo_init_count;
    uint32_t motion_estimation_fifo_init_count;
    uint32_t initial_rate_control_fifo_init_count;
    uint32_t picture_demux_fifo_init_count;
    uint32_t rate_control_tasks_fifo_init_count;
    uint32_t rate_control_fifo_init_count;
    uint32_t mode_decision_configuration_fifo_init_count;
    uint32_t enc_dec_fifo_init_count;
    uint32_t entropy_coding_fifo_init_count;
    uint32_t dlf_fifo_init_count;
    uint32_t cdef_fifo_init_count;
    uint32_t rest_fifo_init_count;

    /*!< Thread count for each process */
    uint32_t picture_analysis_process_init_count;
    uint32_t motion_estimation_process_init_count;
    uint32_t source_based_operations_process_init_count;
    uint32_t mode_decision_configuration_process_init_count;
    uint32_t enc_dec_process_init_count;
    uint32_t entropy_coding_process_init_count;
    uint32_t dlf_process_init_count;
    uint32_t cdef_process_init_count;
    uint32_t rest_process_init_count;
    uint32_t total_process_init_count;

} SequenceControlSet;

typedef struct EbSequenceControlSetInitData {
    EncodeContext *encode_context_ptr;
    int32_t        sb_size;
} EbSequenceControlSetInitData;

typedef struct EbSequenceControlSetInstance {
    EbDctor             dctor;
    EncodeContext *     encode_context_ptr;
    SequenceControlSet *scs_ptr;
    EbHandle            config_mutex;
} EbSequenceControlSetInstance;

/**************************************
     * Extern Function Declarations
     **************************************/
extern EbErrorType eb_sequence_control_set_creator(EbPtr *object_dbl_ptr,
                                                   EbPtr  object_init_data_ptr);

extern EbErrorType eb_sequence_control_set_ctor(SequenceControlSet *object,
                                                EbPtr               object_init_data_ptr);

extern EbErrorType copy_sequence_control_set(SequenceControlSet *dst, SequenceControlSet *src);

extern EbErrorType eb_sequence_control_set_instance_ctor(EbSequenceControlSetInstance *object_ptr);

extern EbErrorType sb_params_init(SequenceControlSet *scs_ptr);

extern EbErrorType derive_input_resolution(EbInputResolution *input_resolution,
                                           uint32_t           input_size);

EbErrorType sb_geom_init(SequenceControlSet *scs_ptr);

#ifdef __cplusplus
}
#endif
#endif // EbSequenceControlSet_h
