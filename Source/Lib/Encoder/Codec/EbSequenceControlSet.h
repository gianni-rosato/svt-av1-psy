// clang-format off
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
    typedef struct SequenceControlSet
    {
        EbDctor                                 dctor;
        EbSvtAv1EncConfiguration                static_config;

        // Encoding Context
        EncodeContext                        *encode_context_ptr;

        // Profile & ID
        uint32_t                                sps_id;
        uint32_t                                vps_id;
        uint32_t                                profile_space;
        uint32_t                                profile_idc;
        uint32_t                                level_idc;
        uint32_t                                tier_idc;
        uint32_t                                chroma_format_idc;
        uint32_t                                max_temporal_layers;
        uint32_t                                bits_for_picture_order_count;
        uint16_t                                subsampling_x;            // add chroma subsampling parameters
        uint16_t                                subsampling_y;

        // Picture deminsions
        uint16_t                                max_input_luma_width;
        uint16_t                                max_input_luma_height;
        uint16_t                                max_input_chroma_width;
        uint16_t                                max_input_chroma_height;
        uint16_t                                max_input_pad_bottom;
        uint16_t                                max_input_pad_right;

        uint32_t                                chroma_width;
        uint32_t                                chroma_height;
        uint32_t                                pad_right;
        uint32_t                                pad_bottom;
        uint16_t                                left_padding;
        uint16_t                                top_padding;
        uint16_t                                right_padding;
        uint16_t                                bot_padding;
        uint32_t                                frame_rate;
        uint32_t                                encoder_bit_depth;

        // Cropping Definitions
        int32_t                                 cropping_left_offset;
        int32_t                                 cropping_right_offset;
        int32_t                                 cropping_top_offset;
        int32_t                                 cropping_bottom_offset;

        // Conformance Window flag
        uint32_t                                conformance_window_flag;

        // Group of Pictures (GOP) Structure
        uint32_t                                max_ref_count;            // Maximum number of reference pictures, however each pred
                                                            //   entry can be less.
        PredictionStructure                  *pred_struct_ptr;
        int32_t                                 intra_period_length;      // The frequency of intra pictures
        uint32_t                                intra_refresh_type;       // 1: CRA, 2: IDR

        // LCU
        uint8_t                                 sb_sz; // lcu_size
        uint8_t                                 max_sb_depth;
        // Coding unit
        uint8_t                                 max_cu_size;
        uint8_t                                 min_cu_size;
        uint8_t                                 max_intra_size;
        uint8_t                                 min_intra_size;
        EbBool                                  intra4x4_flag;
        // Rate Control
        uint32_t                                rate_control_mode;
        uint32_t                                target_bitrate;
        uint32_t                                available_bandwidth;

        // Quantization
       /* uint32_t                                qp;*/
        // Picture Analysis
        uint32_t                                picture_analysis_number_of_regions_per_width;
        uint32_t                                picture_analysis_number_of_regions_per_height;

        // Segments
        uint32_t                                me_segment_column_count_array[MAX_TEMPORAL_LAYERS];
        uint32_t                                me_segment_row_count_array[MAX_TEMPORAL_LAYERS];
        uint32_t                                enc_dec_segment_col_count_array[MAX_TEMPORAL_LAYERS];
        uint32_t                                enc_dec_segment_row_count_array[MAX_TEMPORAL_LAYERS];
        uint32_t                                cdef_segment_column_count;
        uint32_t                                cdef_segment_row_count;

        uint32_t                                rest_segment_column_count;
        uint32_t                                rest_segment_row_count;
        uint32_t                                tf_segment_column_count;
        uint32_t                                tf_segment_row_count;
        EbBool                                  enable_altrefs;
        uint32_t                                scd_delay; //Number of delay frames needed to implement future window for algorithms such as SceneChange or TemporalFiltering
        // Buffers
        uint32_t                                picture_control_set_pool_init_count;
        uint32_t                                picture_control_set_pool_init_count_child;
        uint32_t                                pa_reference_picture_buffer_init_count;
        uint32_t                                reference_picture_buffer_init_count;
        uint32_t                                input_buffer_fifo_init_count;
        uint32_t                                overlay_input_picture_buffer_init_count;
        uint32_t                                output_stream_buffer_fifo_init_count;
        uint32_t                                output_recon_buffer_fifo_init_count;
        uint32_t                                resource_coordination_fifo_init_count;
        uint32_t                                picture_analysis_fifo_init_count;
        uint32_t                                picture_decision_fifo_init_count;
        uint32_t                                motion_estimation_fifo_init_count;
        uint32_t                                initial_rate_control_fifo_init_count;
        uint32_t                                picture_demux_fifo_init_count;
        uint32_t                                rate_control_tasks_fifo_init_count;
        uint32_t                                rate_control_fifo_init_count;
        uint32_t                                mode_decision_configuration_fifo_init_count;
        uint32_t                                enc_dec_fifo_init_count;
        uint32_t                                entropy_coding_fifo_init_count;
        uint32_t                                dlf_fifo_init_count;
        uint32_t                                cdef_fifo_init_count;
        uint32_t                                rest_fifo_init_count;

        uint32_t                                picture_analysis_process_init_count;
        uint32_t                                motion_estimation_process_init_count;
        uint32_t                                source_based_operations_process_init_count;
        uint32_t                                mode_decision_configuration_process_init_count;
        uint32_t                                enc_dec_process_init_count;
        uint32_t                                entropy_coding_process_init_count;
        uint32_t                                dlf_process_init_count;
        uint32_t                                cdef_process_init_count;
        uint32_t                                rest_process_init_count;
        uint32_t                                total_process_init_count;

        uint16_t                                film_grain_random_seed;
        SbParams                               *sb_params_array;
        uint8_t                                 picture_width_in_sb;
        uint8_t                                 picture_height_in_sb;
        uint16_t                                sb_total_count;
        uint16_t                                sb_size_pix;  //sb size in pixels 64/128
        uint16_t                                sb_tot_cnt;   // sb total number
        uint16_t                                max_block_cnt;
        SbGeom                                 *sb_geom;

        EbInputResolution                       input_resolution;
        EbScdMode                               scd_mode;
        EbPmMode                                pm_mode;

        /* MRP (mm-signal; 0: MRP mode 0, 1: MRP mode 1)
        *
        * Default is 0. */
        uint8_t                                 mrp_mode;

        /* CDF (mm-signal; 0: CDF update, 1: no CDF update)
        *
        * Default is 0. */
        uint8_t                                 cdf_mode;

        /* NSQ present (mm-signal; 0: NSQ absent, 1: NSQ present)
        *
        * Default is 1. */
        uint8_t                                 nsq_present;

        /* Down-sampling method @ ME and alt-ref temporal filtering (mm-signal; 0: filtering, 1: decimation)
        *
        * Default is 0. */
        uint8_t                                 down_sampling_method_me_search;
        uint8_t                                 mfmv_enabled; // 1:Enabled  0:Disabled
        uint8_t                                 trans_coeff_shape_array[2][8][4];    // [componantTypeIndex][resolutionIndex][levelIndex][tuSizeIndex]
        EbBlockMeanPrec                         block_mean_calc_prec;
        BitstreamLevel                          level[MAX_NUM_OPERATING_POINTS];
        int32_t                                 film_grain_denoise_strength;
        uint32_t                                reference_count;
        /* over_boundary_block (mm-signal; 0: No over boundary blk allowed, 1: over boundary blk allowed)
        *
        * Default is 0.
        * To enable when md_skip_blk is on*/
        uint8_t                                 over_boundary_block_mode;
        SeqHeader                               seq_header;
        uint8_t                                 compound_mode;
        EbBool                                  use_input_stat_file;
        EbBool                                  use_output_stat_file;
    } SequenceControlSet;

    typedef struct EbSequenceControlSetInitData
    {
        EncodeContext            *encode_context_ptr;
        int32_t                     sb_size;
    } EbSequenceControlSetInitData;

    typedef struct EbSequenceControlSetInstance
    {
        EbDctor                  dctor;
        EncodeContext            *encode_context_ptr;
        SequenceControlSet       *sequence_control_set_ptr;
        EbHandle                    config_mutex;
    } EbSequenceControlSetInstance;

    /**************************************
     * Extern Function Declarations
     **************************************/
    extern EbErrorType eb_sequence_control_set_creator(
        EbPtr *object_dbl_ptr,
        EbPtr object_init_data_ptr);

    extern EbErrorType eb_sequence_control_set_ctor(
        SequenceControlSet* object,
        EbPtr  object_init_data_ptr);

    extern EbErrorType copy_sequence_control_set(
        SequenceControlSet *dst,
        SequenceControlSet *src);

    extern EbErrorType eb_sequence_control_set_instance_ctor(
        EbSequenceControlSetInstance *object_dbl_ptr);

    extern EbErrorType sb_params_init(
        SequenceControlSet *sequence_control_set_ptr);

    extern EbErrorType derive_input_resolution(
        SequenceControlSet *sequence_control_set_ptr,
        uint32_t              input_size);

    EbErrorType sb_geom_init(SequenceControlSet *sequence_control_set_ptr);

#ifdef __cplusplus
}
#endif
#endif // EbSequenceControlSet_h
// clang-format on
