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
#include "EbThreads.h"
#include "EbSystemResourceManager.h"
#include "EbEncodeContext.h"
#include "EbPredictionStructure.h"
#include "noise_model.h"


#ifdef __cplusplus
extern "C" {
#endif
    /************************************
     * Sequence Control Set
     ************************************/
    typedef struct SequenceControlSet_s
    {
        EbSvtAv1EncConfiguration                static_config;

        // Encoding Context
        EncodeContext_t                        *encode_context_ptr;

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
        int32_t                                 subsampling_x;            // add chroma subsampling parameters
        int32_t                                 subsampling_y;

        // Picture deminsions
        uint16_t                                max_input_luma_width;
        uint16_t                                max_input_luma_height;
        uint16_t                                max_input_chroma_width;
        uint16_t                                max_input_chroma_height;
        uint16_t                                max_input_pad_bottom;
        uint16_t                                max_input_pad_right;
        
        uint16_t                                luma_width;
        uint16_t                                luma_height;
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

        // Bitdepth
        EB_BITDEPTH                             input_bitdepth;
        EB_BITDEPTH                             output_bitdepth;

        // Group of Pictures (GOP) Structure
        uint32_t                                max_ref_count;            // Maximum number of reference pictures, however each pred
                                                            //   entry can be less.
        PredictionStructure_t                  *pred_struct_ptr;
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

        uint32_t                                general_progressive_source_flag;
        uint32_t                                general_interlaced_source_flag;
        uint32_t                                general_frame_only_constraint_flag;

        // Rate Control
        uint32_t                                rate_control_mode;
        uint32_t                                target_bitrate;
        uint32_t                                available_bandwidth;

        // Quantization
        uint32_t                                qp;

        // tmvp enable
        uint32_t                                enable_tmvp_sps;

        // MV merge
        uint32_t                                mv_merge_total_count;

        // Picture Analysis
        uint32_t                                picture_analysis_number_of_regions_per_width;
        uint32_t                                picture_analysis_number_of_regions_per_height;

        // Segments
        uint32_t                                me_segment_column_count_array[MAX_TEMPORAL_LAYERS];
        uint32_t                                me_segment_row_count_array[MAX_TEMPORAL_LAYERS];
        uint32_t                                enc_dec_segment_col_count_array[MAX_TEMPORAL_LAYERS];
        uint32_t                                enc_dec_segment_row_count_array[MAX_TEMPORAL_LAYERS];
#if CDEF_M
        uint32_t                                cdef_segment_column_count;
        uint32_t                                cdef_segment_row_count;
#endif
#if REST_M
        uint32_t                                rest_segment_column_count;
        uint32_t                                rest_segment_row_count;
#endif
        // Buffers
        uint32_t                                picture_control_set_pool_init_count;
        uint32_t                                picture_control_set_pool_init_count_child;
        uint32_t                                pa_reference_picture_buffer_init_count;
        uint32_t                                reference_picture_buffer_init_count;
        uint32_t                                input_buffer_fifo_init_count;
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
#if FILT_PROC
        uint32_t                                dlf_fifo_init_count;
        uint32_t                                cdef_fifo_init_count;
        uint32_t                                rest_fifo_init_count;
#endif
        uint32_t                                picture_analysis_process_init_count;
        uint32_t                                motion_estimation_process_init_count;
        uint32_t                                source_based_operations_process_init_count;
        uint32_t                                mode_decision_configuration_process_init_count;
        uint32_t                                enc_dec_process_init_count;
        uint32_t                                entropy_coding_process_init_count;
#if FILT_PROC
        uint32_t                                dlf_process_init_count;
        uint32_t                                cdef_process_init_count;
        uint32_t                                rest_process_init_count;
#endif
        uint32_t                                total_process_init_count;
        
        uint16_t                                film_grain_random_seed;
        SbParams_t                             *sb_params_array;
        uint8_t                                 picture_width_in_sb;
        uint8_t                                 picture_height_in_sb;
        uint16_t                                sb_total_count;
        uint16_t                                sb_size_pix;  //sb size in pixels 64/128
        uint16_t                                sb_tot_cnt;   // sb total number
        uint16_t                                max_block_cnt;
        SbGeom_t                               *sb_geom;

        EbInputResolution                       input_resolution;
        EbScdMode                               scd_mode;
        EbPmMode                                pm_mode;
        uint8_t                                 trans_coeff_shape_array[2][8][4];    // [componantTypeIndex][resolutionIndex][levelIndex][tuSizeIndex]
        EbBlockMeanPrec                         block_mean_calc_prec;

        int32_t                                 num_bits_width;
        int32_t                                 num_bits_height;
        int32_t                                 frame_id_numbers_present_flag;
        int32_t                                 frame_id_length;
        int32_t                                 delta_frame_id_length;
        block_size                               sb_size;                            // Size of the superblock used for this frame
        int32_t                                 mib_size;                           // Size of the superblock in units of MI blocks
        int32_t                                 mib_size_log2;                      // Log 2 of above.
        int32_t                                 order_hint_bits_minus1;
        int32_t                                 force_screen_content_tools;         // 0 - force off
                                                                                    // 1 - force on
                                                                                    // 2 - adaptive
        int32_t                                 force_integer_mv;                   // 0 - Not to force. MV can be in 1/4 or 1/8
                                                                                    // 1 - force to integer
                                                                                    // 2 - adaptive
        int32_t                                 monochrome;
        int32_t                                 enable_filter_intra;                // enables/disables filterintra
        int32_t                                 enable_intra_edge_filter;           // enables/disables corner/edge/upsampling
        int32_t                                 enable_interintra_compound;         // enables/disables interintra_compound
        int32_t                                 enable_masked_compound;             // enables/disables masked compound
        int32_t                                 enable_dual_filter;                 // 0 - disable dual interpolation filter
                                                                                    // 1 - enable vertical and horiz filter selection
        int32_t                                 enable_order_hint;                  // 0 - disable order hint, and related tools:
                                                                                    // jnt_comp, ref_frame_mvs, frame_sign_bias
                                                                                    // if 0, enable_jnt_comp must be set zs 0.
        int32_t                                 enable_jnt_comp;                    // 0 - disable joint compound modes
                                                                                    // 1 - enable it
        int32_t                                 enable_ref_frame_mvs;               // 0 - disable ref frame mvs
                                                                                    // 1 - enable it
        int32_t                                 enable_superres;                    // 0 - Disable superres for the sequence, and disable
                                                                                    //     transmitting per-frame superres enabled flag.
                                                                                    // 1 - Enable superres for the sequence, and also
                                                                                    //     enable per-frame flag to denote if superres is
                                                                                    //     enabled for that frame.
        int32_t                                 enable_cdef;                        // To turn on/off CDEF
        int32_t                                 enable_restoration;                 // To turn on/off loop restoration

        int32_t                                 operating_point_idc[MAX_NUM_OPERATING_POINTS];
        BitstreamLevel                          level[MAX_NUM_OPERATING_POINTS];
        int32_t                                 tier[MAX_NUM_OPERATING_POINTS];
        int32_t                                 reduced_still_picture_hdr;
        int32_t                                 still_picture;
        int32_t                                 timing_info_present;
        int32_t                                 operating_points_decoder_model_cnt;

#if AV1_UPGRADE
        int32_t                                 decoder_model_info_present_flag;
        int32_t                                 display_model_info_present_flag;
#endif
        int32_t                                 film_grain_denoise_strength;
        int32_t                                 film_grain_params_present;  // To turn on/off film grain (on a sequence basis)

#if BASE_LAYER_REF
        int32_t                                 extra_frames_to_ref_islice;
        int32_t                                 max_frame_window_to_ref_islice;
#endif

#if ADP_STATS_PER_LAYER
        uint64_t                                total_count[5];
        uint64_t                                sq_search_count[5];
        uint64_t                                sq_non4_search_count[5];
        uint64_t                                mdc_count[5];
        uint64_t                                pred_count[5];
        uint64_t                                pred1_nfl_count[5];
#endif
    } SequenceControlSet_t;

    typedef struct EbSequenceControlSetInitData_s
    {
        EncodeContext_t            *encode_context_ptr;
        int32_t                     sb_size;
    } EbSequenceControlSetInitData_t;

    typedef struct EbSequenceControlSetInstance_s
    {
        EncodeContext_t            *encode_context_ptr;
        SequenceControlSet_t       *sequence_control_set_ptr;
        EbHandle                    config_mutex;

    } EbSequenceControlSetInstance_t;

    /**************************************
     * Extern Function Declarations
     **************************************/
    extern EbErrorType eb_sequence_control_set_ctor(
        EbPtr *object_dbl_ptr,
        EbPtr  object_init_data_ptr);

    extern EbErrorType copy_sequence_control_set(
        SequenceControlSet_t *dst,
        SequenceControlSet_t *src);

    extern EbErrorType eb_sequence_control_set_instance_ctor(
        EbSequenceControlSetInstance_t **object_dbl_ptr);

    extern EbErrorType sb_params_ctor(
        SequenceControlSet_t *sequence_control_set_ptr);

    extern EbErrorType sb_params_init(
        SequenceControlSet_t *sequence_control_set_ptr);

    extern EbErrorType derive_input_resolution(
        SequenceControlSet_t *sequence_control_set_ptr,
        uint32_t              input_size);

    EbErrorType sb_geom_init(SequenceControlSet_t *sequence_control_set_ptr);

#ifdef __cplusplus
}
#endif
#endif // EbSequenceControlSet_h